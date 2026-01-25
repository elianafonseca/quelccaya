var area = ee.FeatureCollection("projects/ee-egeescripts/assets/quelccaya/quelccaya_4800");
Map.centerObject(area, 12);
var AOI = area;

// --- 1. MODELO DE ELEVAÇÃO (NASADEM) ---
var elevation = ee.Image('NASA/NASADEM_HGT/001').select('elevation');
Map.addLayer(elevation.clip(AOI), {min:4500, max:5500, palette:["green", "lime", "yellow", "orange", "red"]}, 'elevation_AOI');

// --- 2. SELEÇÃO DE IMAGENS (APENAS JUNHO) ---
var imageList = ee.List([
  'COPERNICUS/S2_HARMONIZED/20200618T145731_20200618T150227_T19LBE',
  'COPERNICUS/S2_HARMONIZED/20200618T145731_20200618T150227_T19LCE',
  'COPERNICUS/S2_HARMONIZED/20210623T145731_20210623T150555_T19LBE',
  'COPERNICUS/S2_HARMONIZED/20210623T145731_20210623T150555_T19LCE',
  'COPERNICUS/S2_HARMONIZED/20220628T145741_20220628T150513_T19LBE',
  'COPERNICUS/S2_HARMONIZED/20220628T145741_20220628T150513_T19LCE',
  'COPERNICUS/S2_HARMONIZED/20230628T145729_20230628T145830_T19LBE',
  'COPERNICUS/S2_HARMONIZED/20230628T145729_20230628T145830_T19LCE',
  'COPERNICUS/S2_HARMONIZED/20240607T145731_20240607T150237_T19LBE',
  'COPERNICUS/S2_HARMONIZED/20240607T145731_20240607T150237_T19LCE'
]);

var s2 = ee.ImageCollection('COPERNICUS/S2_HARMONIZED');
var sel = s2.filter(ee.Filter.inList('system:id', imageList));

// --- 3. FUNÇÕES DE PROCESSAMENTO ---
function mosaicByDate(imcol){
  var imlist = imcol.toList(imcol.size());
  var unique_dates = imlist.map(function(im){
    return ee.Image(im).date().format("YYYY-MM-dd");
  }).distinct();
  var mosaic_imlist = unique_dates.map(function(d){
    d = ee.Date(d);
    var im = imcol.filterDate(d, d.advance(1, "day")).mosaic();
    return im.set("system:time_start", d.millis(), "system:id", d.format("YYYY-MM-dd"));
  });
  return ee.ImageCollection(mosaic_imlist);
}

function add_indices(image) {
  var ndvi = image.expression("(nir - red)/(nir+red)", {nir: image.select("B8"), red: image.select("B4")}).rename("NDVI");
  var ndsi = image.expression("(green - swir1)/(green + swir1)", {green: image.select("B3"), swir1: image.select("B11")}).rename("NDSI");
  var csi = image.expression("(nir - swir2)/(nir + swir2)", {nir: image.select("B8"), swir2: image.select("B12")}).rename("CSI");
  var andsi = image.expression("(csi - ndsi)/(csi + ndsi)", {csi: csi, ndsi: ndsi}).rename("ANDSI");
  return image.addBands([ndvi,ndsi,csi,andsi]);
}

var colecao = mosaicByDate(sel).map(add_indices).map(function(img){return img.clip(area)});

// --- 4. CLASSIFICAÇÃO K-MEANS ---
function cluster(image) {
  var entrada = image.select('NDVI', "ANDSI");
  var training = entrada.sample({region: area, scale: 10, numPixels: 5000, seed:0});
  var kmeans = ee.Clusterer.wekaKMeans(10).train(training);
  return image.addBands(entrada.cluster(kmeans));
}

var colecao_clas = colecao.map(cluster);
var listOfImages = colecao_clas.toList(colecao_clas.size());

// --- 5. CONFIGURAÇÃO DA RECLASSIFICAÇÃO ---
var fromList = [0,1,2,3,4,5,6,7,8,9];
var names = ee.List(['Snow','Ice', 'Bare rock', 'Water', 'Biocrust', 'Wetlands']);
var paletteFinal =['FFFFFF', '808080', 'fc8d59', '0000FF', '00FF00', '008000'];
var final_vis = {'min': 1, 'max': 6, 'palette': paletteFinal};

var toLists = {
  '2020': [6,1,3,5,1,3,3,3,2,5],
  '2021': [6,1,3,5,1,3,3,3,2,5],
  '2022': [6,1,3,5,1,3,3,3,2,5],
  '2023': [6,1,3,5,1,3,3,3,2,5],
  '2024': [6,2,3,5,1,3,3,3,3,5]
};

var years = ['2020', '2021', '2022', '2023', '2024'];

// --- 6. LOOP DE PROCESSAMENTO E ESTATÍSTICAS ---
for (var i = 0; i < 5; i++) {
  var img = ee.Image(listOfImages.get(i));
  var year = years[i];
  
  var remapped = img.remap({
    from: fromList,
    to: toLists[year],
    bandName: 'cluster'
  }).rename('class');

  Map.addLayer(remapped, final_vis, year + ' - June');

  // Cálculo de estatísticas usando Reducer.group
  // Banda 0: Valor (Elevation) | Banda 1: Grupo (Class)
  var statsImage = elevation.clip(area).addBands(remapped);
  
  var reducers = ee.Reducer.mean()
    .combine({reducer2: ee.Reducer.min(), sharedInputs: true})
    .combine({reducer2: ee.Reducer.max(), sharedInputs: true});

  var stats = statsImage.reduceRegion({
    reducer: reducers.group({
      groupField: 1,
      groupName: 'class_id',
    }),
    geometry: area,
    scale: 30,
    maxPixels: 1e10
  });

  // Formatação para que o resultado apareça com nomes no Console
  var groups = ee.List(stats.get('groups'));
  var finalStats = groups.map(function(g) {
    var d = ee.Dictionary(g);
    var id = ee.Number(d.get('class_id')).toInt();
    var name = names.get(id.subtract(1));
    return d.set('Nome_Classe', name).select(['Nome_Classe', 'mean', 'min', 'max']);
  });

  print('Elevation (m) - June ' + year, finalStats);
}

// --- 7. LEGENDA ---
var legend = ui.Panel({style: {position: 'bottom-left', padding: '8px 15px'}});
legend.add(ui.Label({value: 'Land Cover Class', style: {fontWeight: 'bold', fontSize: '18px'}}));
var makeRow = function(color, name) {
  var colorBox = ui.Label({style: {backgroundColor: '#' + color, padding: '8px', margin: '0 0 4px 0'}});
  var description = ui.Label({value: name, style: {margin: '0 0 4px 6px'}});
  return ui.Panel({widgets: [colorBox, description], layout: ui.Panel.Layout.Flow('horizontal')});
};
for (var j = 0; j < 6; j++) {
  legend.add(makeRow(paletteFinal[j], names.get(j).getInfo()));
}
Map.add(legend);
