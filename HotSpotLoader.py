import json
import shapefile as shp

def GetHotSpotTrails(FNAME):
    
    shapeobj = shp.Reader(FNAME)
        
    for i in range(len(shapeobj.fields)):
        if 'Age' in shapeobj.fields[i]:
            AgeInd = i-1
        elif 'AgeError' in shapeobj.fields[i]:
            AgeErrorInd = i-1
        elif 'Chain' in shapeobj.fields[i]:
            NameInd = i-1
        elif 'PLATEID1' in shapeobj.fields[i]:
            PlateIDInd = i-1
        
    chains = {}
    
    for shape, record in zip(shapeobj.shapes(), shapeobj.records()):
        
        name = record[NameInd]
        
        # Make the volcanic chain if necessary
        if name not in chains:
            chains[name] = {'lon': [], 'lat': [], 'age': [], 'age_error': [], 'plateid': []}
    
        chains[name]['lon'].append(shape.points[0][0])
        chains[name]['lat'].append(shape.points[0][1])
        chains[name]['age'].append(record[AgeInd])
        chains[name]['age_error'].append(record[AgeErrorInd])
        chains[name]['plateid'].append(record[PlateIDInd])
    
    # Sort the data in terms of age -------
    import itertools as it
    
    for name in chains:
    
        if len(chains[name]['age']) < 2:
            continue
    
        sorted_lists = sorted(it.izip(chains[name]['lon'], chains[name]['lat'],
                              chains[name]['age'], chains[name]['age_error'],
                              chains[name]['plateid']),
                              key=lambda x: x[2])
    
        chains[name]['lon'], chains[name]['lat'], chains[name]['age'], \
            chains[name]['age_error'],chains[name]['plateid'] = it.izip(*sorted_lists)

    return chains


def GetHotSpotLocations(FNAME):
    
    shapeobj = shp.Reader(FNAME)
    
    for i in range(len(shapeobj.fields)):
        if 'NAME' in shapeobj.fields[i]:
            NameInd = i-1

    # Create the dictionary
    hotspots = {}
    for shape, record in zip(shapeobj.shapes(), shapeobj.records()):

        name = record[NameInd]

        hotspots[name] = {}
        hotspots[name]['lon'] = shape.points[0][0]
        hotspots[name]['lat'] = shape.points[0][1]

    return hotspots



def GetHotSpotTrailsFromGeoJSON(FNAME):

    with open(FNAME) as f:
        data = json.load(f)

    chains = {}

    for feature in data['features']:
        name = feature['properties']['Chain']

        # Make the volcanic chain if necessary
        if name not in chains:
            chains[name] = {'lon': [], 'lat': [], 'age': [], 'age_error': [], 'plateid': []}

        chains[name]['lon'].append(feature['geometry']['coordinates'][0])
        chains[name]['lat'].append(feature['geometry']['coordinates'][1])
        chains[name]['age'].append(feature['properties']['Age'])
        chains[name]['age_error'].append(feature['properties']['AgeError'])
        chains[name]['plateid'].append(feature['properties']['PLATEID1'])

    # Sort the data in terms of age -------
    import itertools as it

    for name in chains:

        if len(chains[name]['age']) < 2:
            continue

        sorted_lists = sorted(it.izip(chains[name]['lon'], chains[name]['lat'],
                              chains[name]['age'], chains[name]['age_error'],
                              chains[name]['plateid']),
                              key=lambda x: x[2])

        chains[name]['lon'], chains[name]['lat'], chains[name]['age'], \
            chains[name]['age_error'],chains[name]['plateid'] = it.izip(*sorted_lists)

    return chains

def GetHotSpotLocationsFromGeoJSON(FNAME):

    with open(FNAME) as f:
        data = json.load(f)

    # Create the dictionary
    hotspots = {}
    for feature in data['features']:

        name = feature['properties']['NAME']

        hotspots[name] = {}
        hotspots[name]['lon'] = feature['geometry']['coordinates'][0]
        hotspots[name]['lat'] = feature['geometry']['coordinates'][1]

    return hotspots

