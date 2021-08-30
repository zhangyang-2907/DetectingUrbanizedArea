//----------------------------------------------------------------------------------------
//                      Step 1: calculating annual mean BCI
//----------------------------------------------------------------------------------------
//Note: there are two places that need to set the temporal periods manually. 
//One is in Data section (Lines 89-97) and the other is in Annual mean BCI calculation (Line 232).

//________________________________________________________________________________________
//                                study area (Taken Fuzhou city as an example)
//________________________________________________________________________________________
var roi = /* color: #fcfafa */ee.Geometry.Polygon(
        [[[119.09, 26.18],
          [119.53, 26.18],
          [119.53, 25.845],
          [119.09, 25.845]]]);  
          
Map.centerObject(roi, 10);

//________________________________________________________________________________________
//                                      functions
//________________________________________________________________________________________
//(1) Masking cloud pixels in Landsat SR data
var maskClouds_SR = function(image) {
                var qa = ee.Image(image).select('pixel_qa');   
                var mask = qa.bitwiseAnd(2).eq(2).or(qa.bitwiseAnd(4).eq(4));
                var clean = ee.Image(image).updateMask(mask);
                return clean;
                };

//(2) Using the masks created from the Landsat SR data to mask the cloud pixels in Landsat TOA or raw data
var dataPreprocessing = function (data_SR,data_TOA,startDate,endDate){
  
    var landsat_SR = ee.ImageCollection(data_SR).filterBounds(roi).filterDate(startDate,endDate);
    var landsat_TOA = ee.ImageCollection(data_TOA).filterBounds(roi).filterDate(startDate,endDate);
    
    // Define an inner join.
    var innerJoin = ee.Join.inner();

    // Specify an equals filter for image timestamps.
    var filterTimeEq = ee.Filter.equals({leftField: 'system:time_start',rightField: 'system:time_start'});

    // Apply the join.
    var innerJoined_landsat = innerJoin.apply(landsat_SR, landsat_TOA, filterTimeEq);
    
    // Map a function to merge the results in the output FeatureCollection.
    var joined_Landsat = innerJoined_landsat.map(function(feature) {
        return ee.Image.cat(feature.get('primary'), feature.get('secondary'));
    });

    //masking the clouds
    var landsat_imageNoClouds = ee.ImageCollection(joined_Landsat.map(maskClouds_SR));
    
    //output the TOA bands used for TC transformation
    return landsat_imageNoClouds.map(function(img){return img.clip(roi);});
};

//(3) Calculating the three TC components
var calculateTasseledCap = function (image,b) {
    var TC1 = image.expression('(B * BRIGHTNESS)',{'B':b, 'BRIGHTNESS': brightness_coefficents});
    var TC2 = image.expression('(B * GREENNESS)',{'B':b,'GREENNESS': greenness_coefficents});
    var TC3 = image.expression('(B * WETNESS)',{'B':b,'WETNESS': wetness_coefficents});
    
    TC1 = TC1.reduce(ee.call("Reducer.sum")).rename('TC1');
    TC2 = TC2.reduce(ee.call("Reducer.sum")).rename('TC2');
    TC3 = TC3.reduce(ee.call("Reducer.sum")).rename('TC3');
	
    var tasseled_cap = ee.Image(TC1).addBands(TC2).addBands(TC3)
                       .set('system:time_start',image.get('system:time_start'));
    return tasseled_cap;
};

//(4) Calculating the annual means of BCI indicator
var annualMean = function (dataSet,indicator,startYear,endYear,startMonth,durationYear){
    var years = ee.List.sequence(startYear,endYear);  
    var meanAnnual_indi = years.map(function(y) {
        var start = ee.Date.fromYMD(y, startMonth, 1);  
        var stop = start.advance(durationYear, 'year');  
        var collection = dataSet.select([indicator]).filterDate(start, stop); 
        var mean = collection.reduce(ee.Reducer.mean()).set('year', y);
        return mean.addBands(mean.metadata('year'));
    });
    return ee.ImageCollection(meanAnnual_indi);
    
};

//________________________________________________________________________________________
//                  data (******need to set the temporal period manually******)
//________________________________________________________________________________________

var l8Col_TOA_masked = dataPreprocessing("LANDSAT/LC08/C01/T1_SR","LANDSAT/LC08/C01/T1_TOA","2010-1-1","2019-12-31");
//print(l8Col_TOA_masked);
//Map.addLayer(l8Col_TOA_masked.first());

var l7Col_TOA_masked = dataPreprocessing("LANDSAT/LE07/C01/T1_SR","LANDSAT/LE07/C01/T1_TOA","2010-1-1","2019-12-31");
//print(l7Col_TOA_masked);
//Map.addLayer(l7Col_TOA_masked);

var l5Col_DN_masked = dataPreprocessing("LANDSAT/LT05/C01/T1_SR","LANDSAT/LT05/C01/T1","2010-1-1","2019-12-31");
//print(l5Col_DN_masked);
//Map.addLayer(l5Col_DN_masked);

//________________________________________________________________________________________
//                        TC calculation for Landsat 8 TOA
//________________________________________________________________________________________

//Coefficients are only for Landsat 8 TOA
var brightness_coefficents = ee.Image([0.3029, 0.2786, 0.4733, 0.5599, 0.508, 0.1872]);
var greenness_coefficents = ee.Image([-0.2941, -0.243, -0.5424, 0.7276, 0.0713, -0.1608]);
var wetness_coefficents = ee.Image([0.1511, 0.1973, 0.3283, 0.3407, -0.7117, -0.4559]);

var l8_tasseled_cap = l8Col_TOA_masked.map(function(image) {
  //masking water pixels first using AWEI
  var awei_sh = image.expression('((BLUE+2.5*GREEN-1.5*(NIR+SWIR1)-0.25*SWIR2)/10000)',
              {
               'NIR': image.select('B5'),'BLUE': image.select('B2'),'SWIR1': image.select('B6'),
               'SWIR2': image.select('B7'),'GREEN': image.select('B3')
              });
              
  var water_mask =  awei_sh.lt(0);           
  var image_waterMasked = image.updateMask(water_mask);
  
  var b = image_waterMasked.select('B2_1','B3_1','B4_1','B5_1','B6_1','B7_1');//TOA
  var l8_tc = calculateTasseledCap(image_waterMasked, b);
  return l8_tc;
});
                          
//print(l8_tasseled_cap);
//Map.addLayer(l8_tasseled_cap,{},'Landsat 8 Tasseled Cap');

//________________________________________________________________________________________
//                        TC calculation for Landsat 7 TOA
//________________________________________________________________________________________

//Coefficients are only for Landsat 7 TOA
brightness_coefficents= ee.Image([0.3561, 0.3972, 0.3904, 0.6966, 0.2286, 0.1596]);
greenness_coefficents= ee.Image([-0.3344, -0.3544, -0.4556, 0.6966, -0.0242, -0.2630]);
wetness_coefficents= ee.Image([0.2626, 0.2141, 0.0926, 0.0656, -0.7629, -0.5388]); 

var l7_tasseled_cap = l7Col_TOA_masked.map(function(image) {
  //masking water pixels first using AWEI
  var awei_sh = image.expression('((BLUE+2.5*GREEN-1.5*(NIR+SWIR1)-0.25*SWIR2)/10000)',
              {
               'NIR': image.select('B4'),'BLUE': image.select('B1'),'SWIR1': image.select('B5'),
               'SWIR2': image.select('B7'),'GREEN': image.select('B2')
              });
              
  var water_mask =  awei_sh.lt(0);           
  var image_waterMasked = image.updateMask(water_mask);
  
  var b = image_waterMasked.select('B1_1','B2_1','B3_1','B4_1','B5_1','B7_1');//TOA
  var l7_tc = calculateTasseledCap(image_waterMasked, b);
  return l7_tc;
});
                          
//print(l7_tasseled_cap);
//Map.addLayer(l7_tasseled_cap,{},'Landsat 7 Tasseled Cap');

//________________________________________________________________________________________
//                        TC calculation for Landsat 5 raw DN
//________________________________________________________________________________________
 
  //Coefficients are only for Landsat 5 raw DN
brightness_coefficents = ee.Image([0.2909, 0.2493, 0.4806, 0.5568, 0.4438, 0.1706]);
greenness_coefficents = ee.Image([-0.2728, -0.2174, -0.5508, 0.7221, 0.0733, -0.1648]);
wetness_coefficents= ee.Image([0.1446, 0.1761, 0.3322, 0.3396, -0.6210, -0.4186]);
  
var l5_tasseled_cap = l5Col_DN_masked.map(function(image) {
  //masking water pixels first using AWEI
  var awei_sh = image.expression('((BLUE+2.5*GREEN-1.5*(NIR+SWIR1)-0.25*SWIR2)/10000)',
              {
               'NIR': image.select('B4'),'BLUE': image.select('B1'),'SWIR1': image.select('B5'),
               'SWIR2': image.select('B7'),'GREEN': image.select('B2')
              });
              
  var water_mask =  awei_sh.lt(0);           
  var image_waterMasked = image.updateMask(water_mask);
  
  var b = image_waterMasked.select('B1_1','B2_1','B3_1','B4_1','B5_1','B7_1');//raw DN 
  var l5_tc = calculateTasseledCap(image_waterMasked, b);
  return l5_tc;
});
                          
//print(l5_tasseled_cap);
//Map.addLayer(l5_tasseled_cap,{},'Landsat 5 Tasseled Cap');  
  
//________________________________________________________________________________________
//                        BCI calculation 
//________________________________________________________________________________________

var L578_TC_merge = ee.ImageCollection(l8_tasseled_cap.merge(l7_tasseled_cap.merge(l5_tasseled_cap)));
//print(L578_TC_merge);
//Map.addLayer(L578_TC_merge);

//Calculating the max and min values of TC bands in each image
var TC_dataSet = L578_TC_merge.map(function(image){
  
  var minMax_TC = image.reduceRegion({reducer:ee.Reducer.minMax(), geometry: roi, maxPixels:1e13, tileScale: 16});
  image = image.set('TC1_max',minMax_TC.get('TC1_max'),'TC1_min',minMax_TC.get('TC1_min'),
                    'TC2_max',minMax_TC.get('TC2_max'),'TC2_min',minMax_TC.get('TC2_min'),
                    'TC3_max',minMax_TC.get('TC3_max'),'TC3_min',minMax_TC.get('TC3_min'));
  return image;
}); 

//print(TC_dataSet);

//Filtering the images without valid pixels
var TC_dataSet_filtered = TC_dataSet.filterMetadata('TC1_max',"not_equals",null);
//print(TC_dataSet_filtered);

//Calculating the BCI value for each image
var BCI_dataSet = TC_dataSet_filtered.map(function(image){
  
  var minMax_TC1 = ee.Number(image.get('TC1_max')).subtract(ee.Number(image.get('TC1_min')));
  var minMax_TC2 = ee.Number(image.get('TC2_max')).subtract(ee.Number(image.get('TC2_min')));
  var minMax_TC3 = ee.Number(image.get('TC3_max')).subtract(ee.Number(image.get('TC3_min')));

  var H = image.select('TC1').subtract(ee.Number(image.get('TC1_min'))).divide(minMax_TC1);
  var V = image.select('TC2').subtract(ee.Number(image.get('TC2_min'))).divide(minMax_TC2);
  var L = image.select('TC3').subtract(ee.Number(image.get('TC3_min'))).divide(minMax_TC3);

  var BCI =  H.add(L).divide(2).subtract(V).divide(H.add(L).divide(2).add(V)).rename('BCI')
             .set('system:time_start',image.get('system:time_start'));
  
  return BCI;
});

print(BCI_dataSet);
//Map.addLayer(BCI_dataSet);

//________________________________________________________________________________________
//   Annual mean BCI calculation (******need to set the temporal period manually******)
//________________________________________________________________________________________
var annualMean_BCI = annualMean(BCI_dataSet,'BCI',2010,2019,1,1);
print(annualMean_BCI);
//Map.addLayer(annualMean_BCI.select('BCI_mean'));

//________________________________________________________________________________________
//        Converting annual mean BCI collection to one image with multiple bands
//________________________________________________________________________________________
var first = ee.Image(0);
//print(first);

//converting collection to image with multiple bands
var convertColtoImg = function(image,previous){
  //var previous = ee.Image(previous);
  var added = ee.Image(previous).addBands(image);
  return ee.Image([added]);
};

var am_BCI_image = annualMean_BCI.iterate(convertColtoImg, first);
print(am_BCI_image);
Map.addLayer(ee.Image(am_BCI_image));

//________________________________________________________________________________________

Export.image.toAsset({
  image: ee.Image(am_BCI_image),
  description: 'BCI_FZ_1019',
  scale: 30,  
  region: roi
});
  






  