// ============================================================
// BMTF manuscript code: GEE data preparation (Section 2)
//
// Export per 10-day window:
//   - TCG.tif
//   - valid.tif
//   - nValid.tif
//   - S1I.tif
//   - RGB.tif
//   - steps_metadata.csv
//
// Notes:
// - This script is clarity-oriented for open-source sharing.
// - Replace asset IDs and export folder before running.
// ============================================================


// ------------------------------ 0. User settings ------------------------------
var ROI_CALC_ID = 'projects/your-project/assets/roi_calc';
var ROI_OUT_ID  = 'projects/your-project/assets/roi_out';

// Replace both placeholders before running:
//   START_DATE_STR = 'YYYY-MM-DD'  (inclusive)
//   END_DATE_STR   = 'YYYY-MM-DD'  (exclusive)
var START_DATE_STR = 'YYYY-MM-DD';   // e.g. '2025-04-01'
var END_DATE_STR   = 'YYYY-MM-DD';   // e.g. '2025-12-01'
var START = ee.Date(START_DATE_STR); // inclusive
var END   = ee.Date(END_DATE_STR);   // exclusive
var STEP_DAYS = 10;

var SCALE = 10;
var MAX_PIXELS = 1e13;
var EXPORT_FOLDER = 'BMTF_GEE_EXPORT';

// S2 quality
var AOT_MAX = 0.6;
var TCG_MIN = -0.12;
var TCG_MAX =  0.35;

// S1 fallback and Lee filter
var S1_SEARCH_DAYS = 25;
var LEE_KERNEL_RADIUS = 1;  // 3x3 window
var LEE_ENL = 4;


// ------------------------------ 1. Collections ------------------------------
var roiCalc = ee.FeatureCollection(ROI_CALC_ID);
var roiOut  = ee.FeatureCollection(ROI_OUT_ID);
var regionGeom = roiOut.geometry();

var S2_SR  = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED');
var S2_L1C = ee.ImageCollection('COPERNICUS/S2_HARMONIZED');
var S1_GRD = ee.ImageCollection('COPERNICUS/S1_GRD');

// Use Sentinel-2 native 10 m projection as reference grid.
var REF_PROJ = ee.Image(
  S2_SR.filterBounds(roiCalc).filterDate(START, END).first()
).select('B2').projection();

function lockToRef(img) {
  return ee.Image(img).reproject(REF_PROJ.atScale(SCALE));
}

var projInfo = REF_PROJ.atScale(SCALE).getInfo();
var REF_CRS = projInfo.crs;
var REF_TRANSFORM = projInfo.transform;
print('REF_CRS:', REF_CRS);
print('REF_TRANSFORM:', REF_TRANSFORM);


// ------------------------------ 2. Sentinel-2 utilities ------------------------------
function qualityMaskL2A(img) {
  img = ee.Image(img);

  var scl = img.select('SCL');
  var sclOk = scl.neq(0)
    .and(scl.neq(1))
    .and(scl.neq(3))
    .and(scl.neq(8))
    .and(scl.neq(9))
    .and(scl.neq(10))
    .and(scl.neq(11));

  var aot = img.select('AOT').multiply(0.001);
  var aotOk = aot.lt(AOT_MAX).unmask(1);

  var key = img.select(['B8', 'B11']).multiply(0.0001);
  var keyOk = key.unmask(0).gt(0).reduce(ee.Reducer.min());

  return sclOk.and(aotOk).and(keyOk);
}

function tcgFromL1C(l1c) {
  l1c = ee.Image(l1c);
  var r = l1c.select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12']).multiply(0.0001).toFloat();

  var B2 = r.select('B2');
  var B3 = r.select('B3');
  var B4 = r.select('B4');
  var B8 = r.select('B8');
  var B11 = r.select('B11');
  var B12 = r.select('B12');

  return r.expression(
    '-0.2848*B2 - 0.2435*B3 - 0.5436*B4 + 0.7243*B8 + 0.0840*B11 - 0.1800*B12',
    {B2: B2, B3: B3, B4: B4, B8: B8, B11: B11, B12: B12}
  ).rename('TCG');
}

function buildS2Step(startDate, endDate, labelStr, idxNum) {
  var srCol  = S2_SR.filterBounds(roiCalc).filterDate(startDate, endDate);
  var l1cCol = S2_L1C.filterBounds(roiCalc).filterDate(startDate, endDate);

  var join = ee.Join.saveFirst('l1c');
  var filt = ee.Filter.equals({leftField: 'system:index', rightField: 'system:index'});
  var srJoined = ee.ImageCollection(join.apply(srCol, l1cCol, filt));
  var paired = srJoined.filter(ee.Filter.notNull(['l1c']));
  var nImgPair = paired.size();

  var midMs = startDate.millis().add(endDate.millis()).divide(2);
  var s2CenterMs = ee.Number(ee.Algorithms.If(
    nImgPair.gt(0),
    ee.Number(paired.aggregate_mean('system:time_start')),
    midMs
  ));

  var perImg = paired.map(function(srImg) {
    srImg = ee.Image(srImg);
    var l1c = ee.Image(srImg.get('l1c'));
    var valid = qualityMaskL2A(srImg).rename('valid');
    var tcg = tcgFromL1C(l1c).updateMask(valid).clamp(TCG_MIN, TCG_MAX).rename('TCG');

    return ee.Image.cat([
      tcg,
      valid.toUint8().rename('valid')
    ]).copyProperties(srImg, ['system:time_start']);
  });

  var tcgMed = ee.Image(ee.Algorithms.If(
    nImgPair.gt(0),
    perImg.select('TCG').median().rename('TCG'),
    ee.Image(0).rename('TCG').updateMask(ee.Image(0))
  ));

  var nValid = ee.Image(ee.Algorithms.If(
    nImgPair.gt(0),
    perImg.select('valid').sum().rename('nValid').toUint16().unmask(0),
    ee.Image(0).rename('nValid').toUint16().unmask(0)
  ));

  var valid = nValid.gte(1).rename('valid').toUint8();

  var out = ee.Image.cat([
    lockToRef(tcgMed).rename('TCG'),
    lockToRef(valid).rename('valid'),
    lockToRef(nValid).rename('nValid')
  ]).clip(roiOut);

  return out.set({
    idx: idxNum,
    label: labelStr,
    nImgPair: nImgPair,
    s2_center_ms: s2CenterMs
  });
}

function buildRGBStep(startDate, endDate) {
  var srCol = S2_SR.filterBounds(roiCalc).filterDate(startDate, endDate);
  var rgbCol = srCol.map(function(im) {
    im = ee.Image(im);
    var valid = qualityMaskL2A(im);
    return im.select(['B4', 'B3', 'B2']).multiply(0.0001).toFloat().updateMask(valid);
  });

  var rgb = ee.Image(ee.Algorithms.If(
    rgbCol.size().gt(0),
    rgbCol.median(),
    ee.Image.constant([0, 0, 0]).rename(['B4', 'B3', 'B2']).updateMask(ee.Image(0))
  ));
  return lockToRef(rgb).clip(roiOut);
}


// ------------------------------ 3. Sentinel-1 utilities ------------------------------
function dbToLinear(img) {
  img = ee.Image(img);
  return ee.Image(10).pow(img.divide(10));
}

function linearToDb(img) {
  img = ee.Image(img);
  return img.max(1e-10).log10().multiply(10);
}

function leeFilterLinear(img, bandName, radius, enl) {
  img = ee.Image(img);
  var band = img.select(bandName);
  var kernel = ee.Kernel.square({radius: radius, units: 'pixels', normalize: false});

  var mean = band.reduceNeighborhood({reducer: ee.Reducer.mean(), kernel: kernel});
  var variance = band.reduceNeighborhood({reducer: ee.Reducer.variance(), kernel: kernel});

  var noiseVar = mean.pow(2).divide(enl);
  var weight = variance.subtract(noiseVar).divide(variance).clamp(0, 1);
  return mean.add(weight.multiply(band.subtract(mean))).rename(bandName);
}

function preprocessS1(img) {
  img = ee.Image(img);
  var vvDb = img.select('VV').rename('VV');
  var vhDb = img.select('VH').rename('VH');

  var vvLin = dbToLinear(vvDb).rename('VV');
  var vhLin = dbToLinear(vhDb).rename('VH');

  var vvLeeLin = leeFilterLinear(vvLin, 'VV', LEE_KERNEL_RADIUS, LEE_ENL);
  var vhLeeLin = leeFilterLinear(vhLin, 'VH', LEE_KERNEL_RADIUS, LEE_ENL);

  var vvLeeDb = linearToDb(vvLeeLin).rename('VV');
  var vhLeeDb = linearToDb(vhLeeLin).rename('VH');

  return ee.Image.cat([vvLeeDb, vhLeeDb]).copyProperties(img, ['system:time_start']);
}

function getS1InStepOrFallback(startDate, endDate, s2CenterMs) {
  var s1win = S1_GRD
    .filterBounds(roiCalc)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.eq('resolution_meters', 10))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .map(preprocessS1);

  var nS1 = s1win.size();
  var midMs = startDate.millis().add(endDate.millis()).divide(2);
  var s1CenterWinMs = ee.Number(ee.Algorithms.If(
    nS1.gt(0),
    ee.Number(s1win.aggregate_mean('system:time_start')),
    midMs
  ));

  var s1iWin = ee.Image(ee.Algorithms.If(
    nS1.gt(0),
    s1win.select('VV').median().add(s1win.select('VH').median()).rename('S1I'),
    ee.Image(0).rename('S1I').updateMask(ee.Image(0))
  ));
  s1iWin = lockToRef(s1iWin).clip(roiOut);

  // Fallback when no SAR in current 10-day window.
  var t0 = ee.Date(s2CenterMs);
  var s1all = S1_GRD
    .filterBounds(roiCalc)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.eq('resolution_meters', 10))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .map(preprocessS1);

  var colBefore = s1all.filterDate(t0.advance(-S1_SEARCH_DAYS, 'day'), t0);
  var colAfter  = s1all.filterDate(t0, t0.advance(S1_SEARCH_DAYS, 'day'));

  var hasBefore = ee.Number(colBefore.size()).gt(0);
  var hasAfter  = ee.Number(colAfter.size()).gt(0);
  var before = colBefore.sort('system:time_start', false).first();
  var after  = colAfter.sort('system:time_start', true).first();

  var big = ee.Number(9e9);
  var dtB = ee.Number(ee.Algorithms.If(
    hasBefore,
    t0.difference(ee.Date(ee.Image(before).get('system:time_start')), 'day').abs(),
    big
  ));
  var dtA = ee.Number(ee.Algorithms.If(
    hasAfter,
    t0.difference(ee.Date(ee.Image(after).get('system:time_start')), 'day').abs(),
    big
  ));

  var useBefore = dtB.lte(dtA);
  var hasNearest = hasBefore.or(hasAfter);
  var nearest = ee.Image(ee.Algorithms.If(useBefore, before, after));

  var s1CenterFbMs = ee.Number(ee.Algorithms.If(
    hasNearest,
    ee.Number(nearest.get('system:time_start')),
    midMs
  ));

  var s1iFb = ee.Image(ee.Algorithms.If(
    hasNearest,
    nearest.select('VV').add(nearest.select('VH')).rename('S1I'),
    ee.Image(0).rename('S1I').updateMask(ee.Image(0))
  ));
  s1iFb = lockToRef(s1iFb).clip(roiOut);

  var fallbackUsed = nS1.eq(0).and(hasNearest);

  return ee.Dictionary({
    img: ee.Image(ee.Algorithms.If(nS1.gt(0), s1iWin, s1iFb)),
    s1_center_ms: ee.Number(ee.Algorithms.If(nS1.gt(0), s1CenterWinMs, s1CenterFbMs)),
    nS1: nS1,
    fallback_used: ee.Number(fallbackUsed)
  });
}


// ------------------------------ 4. Build step list ------------------------------
var nSteps = END.difference(START, 'day').divide(STEP_DAYS).floor();
var stepIdx = ee.List.sequence(0, nSteps.subtract(1));

function stepDict(i) {
  i = ee.Number(i);
  var s = START.advance(i.multiply(STEP_DAYS), 'day');
  var e = s.advance(STEP_DAYS, 'day');
  var label = s.format('YYYYMMdd').cat('_').cat(e.advance(-1, 'day').format('YYYYMMdd'));
  return ee.Dictionary({idx: i, start: s, end: e, label: label});
}

var stepList = stepIdx.map(stepDict);


// ------------------------------ 5. Export helpers ------------------------------
function exportOneStep(iClient) {
  var d = ee.Dictionary(stepList.get(iClient));
  var idx = ee.Number(d.get('idx'));
  var s = ee.Date(d.get('start'));
  var e = ee.Date(d.get('end'));
  var label = ee.String(d.get('label'));

  var s2step = buildS2Step(s, e, label, idx);
  var rgb = buildRGBStep(s, e);
  var s1pack = getS1InStepOrFallback(s, e, ee.Number(s2step.get('s2_center_ms')));
  var s1i = ee.Image(s1pack.get('img'));

  var prefix = ee.String('STEP_').cat(idx.format('%02d')).cat('_').cat(label).getInfo();

  Export.image.toDrive({
    image: s2step.select('TCG').toFloat(),
    description: prefix + '_TCG',
    folder: EXPORT_FOLDER,
    fileNamePrefix: prefix + '/TCG',
    region: regionGeom,
    scale: SCALE,
    crs: REF_CRS,
    crsTransform: REF_TRANSFORM,
    maxPixels: MAX_PIXELS,
    fileFormat: 'GeoTIFF'
  });

  Export.image.toDrive({
    image: s2step.select('valid').toUint8(),
    description: prefix + '_valid',
    folder: EXPORT_FOLDER,
    fileNamePrefix: prefix + '/valid',
    region: regionGeom,
    scale: SCALE,
    crs: REF_CRS,
    crsTransform: REF_TRANSFORM,
    maxPixels: MAX_PIXELS,
    fileFormat: 'GeoTIFF'
  });

  Export.image.toDrive({
    image: s2step.select('nValid').toUint16(),
    description: prefix + '_nValid',
    folder: EXPORT_FOLDER,
    fileNamePrefix: prefix + '/nValid',
    region: regionGeom,
    scale: SCALE,
    crs: REF_CRS,
    crsTransform: REF_TRANSFORM,
    maxPixels: MAX_PIXELS,
    fileFormat: 'GeoTIFF'
  });

  Export.image.toDrive({
    image: s1i.select('S1I').toFloat(),
    description: prefix + '_S1I',
    folder: EXPORT_FOLDER,
    fileNamePrefix: prefix + '/S1I',
    region: regionGeom,
    scale: SCALE,
    crs: REF_CRS,
    crsTransform: REF_TRANSFORM,
    maxPixels: MAX_PIXELS,
    fileFormat: 'GeoTIFF'
  });

  Export.image.toDrive({
    image: rgb.select(['B4', 'B3', 'B2']).toFloat(),
    description: prefix + '_RGB',
    folder: EXPORT_FOLDER,
    fileNamePrefix: prefix + '/RGB',
    region: regionGeom,
    scale: SCALE,
    crs: REF_CRS,
    crsTransform: REF_TRANSFORM,
    maxPixels: MAX_PIXELS,
    fileFormat: 'GeoTIFF'
  });

  print('Created tasks for', prefix);
}

function buildStepMetadata(d) {
  d = ee.Dictionary(d);
  var idx = ee.Number(d.get('idx'));
  var s = ee.Date(d.get('start'));
  var e = ee.Date(d.get('end'));
  var label = ee.String(d.get('label'));

  var s2step = buildS2Step(s, e, label, idx);
  var s1pack = getS1InStepOrFallback(s, e, ee.Number(s2step.get('s2_center_ms')));

  return ee.Feature(null, {
    idx: idx,
    label: label,
    start: s.format('YYYY-MM-dd'),
    end: e.advance(-1, 'day').format('YYYY-MM-dd'),
    nImgPair: ee.Number(s2step.get('nImgPair')),
    s2_center_ms: ee.Number(s2step.get('s2_center_ms')),
    nS1: ee.Number(s1pack.get('nS1')),
    s1_center_ms: ee.Number(s1pack.get('s1_center_ms')),
    fallback_used: ee.Number(s1pack.get('fallback_used'))
  });
}

function exportAll() {
  var n = nSteps.getInfo();
  for (var i = 0; i < n; i++) {
    exportOneStep(i);
  }

  var stepMetaFc = ee.FeatureCollection(stepList.map(buildStepMetadata));
  Export.table.toDrive({
    collection: stepMetaFc,
    description: 'steps_metadata_csv',
    folder: EXPORT_FOLDER,
    fileNamePrefix: 'steps_metadata',
    fileFormat: 'CSV'
  });

  print('All export tasks created. Run them manually in the Tasks panel.');
}


// ------------------------------ 6. Quick preview ------------------------------
function previewStep(iClient) {
  var d = ee.Dictionary(stepList.get(iClient));
  var idx = ee.Number(d.get('idx'));
  var s = ee.Date(d.get('start'));
  var e = ee.Date(d.get('end'));
  var label = ee.String(d.get('label'));

  var s2step = buildS2Step(s, e, label, idx);
  var rgb = buildRGBStep(s, e);
  var s1pack = getS1InStepOrFallback(s, e, ee.Number(s2step.get('s2_center_ms')));
  var s1i = ee.Image(s1pack.get('img'));

  Map.layers().reset([]);
  Map.centerObject(roiOut, 11);
  Map.addLayer(rgb, {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.25}, 'RGB');
  Map.addLayer(s2step.select('TCG'), {min: TCG_MIN, max: TCG_MAX, palette: ['0000ff', '00ff00', 'ff0000']}, 'TCG');
  Map.addLayer(s2step.select('valid'), {min: 0, max: 1, palette: ['000000', '00ff00']}, 'valid', false);
  Map.addLayer(s2step.select('nValid'), {min: 0, max: 5, palette: ['000000', 'ffff00', 'ff0000']}, 'nValid', false);
  Map.addLayer(s1i, {min: -40, max: 0}, 'S1I', false);
  print('Preview:', idx, label);
}


// Preview one step, then call exportAll().
previewStep(0);
// exportAll();
