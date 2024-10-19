__author__ = 'Kozo Nishida'
__email__ = 'kozo.nishida@gmail.com'
__version__ = '0.0.1'
__license__ = 'MIT'

import re
from datetime import datetime

class RmbSpectrum2:
    def __init__(self):
        self.info = {}
        self.data = None  # This will hold the PKPeak DataFrame
        self.collisionEnergy = None
        self.precursorMz = None
        self.rt = None

def addProperty(sp, property_name, property_type):
    # This function seems to be a placeholder for adding properties to the RmbSpectrum2 object
    # In Python, we don't need to pre-define properties, so this function can be omitted
    return sp

def setData(sp, data):
    sp.data = data
    return sp

def parseMbRecord(filename, readAnnotation=True):
    with open(filename, 'r') as fileConnection:
        record = fileConnection.readlines()

    recData = {}
    recData['ACCESSION'] = re.findall(r'ACCESSION:(.+)', ''.join(record))[0].strip()
    recData['RECORD_TITLE'] = re.findall(r'RECORD_TITLE:(.+)', ''.join(record))[0].strip()
    date_str = re.findall(r'DATE:(.+)', ''.join(record))[0].strip()
    recData['DATE'] = datetime.strptime(date_str, '%Y.%m.%d').strftime('%Y.%m.%d')
    recData['AUTHORS'] = re.findall(r'AUTHORS:(.+)', ''.join(record))[0].strip()
    recData['LICENSE'] = re.findall(r'LICENSE:(.+)', ''.join(record))[0].strip()
    recData['COPYRIGHT'] = re.findall(r'COPYRIGHT:(.+)', ''.join(record))[0].strip()

    # COMMENTS
    recData['COMMENT'] = [comment.strip() for comment in re.findall(r'COMMENT:(.+)', ''.join(record))]

    # CH$ fields
    recData['CH$NAME'] = [name.strip() for name in re.findall(r'CH\$NAME:(.+)', ''.join(record))]
    recData['CH$COMPOUND_CLASS'] = re.findall(r'CH\$COMPOUND_CLASS:(.+)', ''.join(record))[0].strip()
    recData['CH$FORMULA'] = re.findall(r'CH\$FORMULA:(.+)', ''.join(record))[0].strip()
    recData['CH$EXACT_MASS'] = float(re.findall(r'CH\$EXACT_MASS:(.+)', ''.join(record))[0].strip())
    recData['CH$SMILES'] = re.findall(r'CH\$SMILES:(.+)', ''.join(record))[0].strip()
    recData['CH$IUPAC'] = re.findall(r'CH\$IUPAC:(.+)', ''.join(record))[0].strip()

    # CH$LINK
    recData['CH$LINK'] = [link.strip() for link in re.findall(r'CH\$LINK:(.+)', ''.join(record))]

    # AC$INSTRUMENT
    recData['AC$INSTRUMENT'] = re.findall(r'AC\$INSTRUMENT:(.+)', ''.join(record))[0].strip()
    recData['AC$INSTRUMENT_TYPE'] = re.findall(r'AC\$INSTRUMENT_TYPE:(.+)', ''.join(record))[0].strip()

    # Version handling (simplified)
    Version = 2
    ac_ms = {}
    ac_ms['MS_TYPE'] = re.findall(r'AC\$MASS_SPECTROMETRY: MS_TYPE (.+)', ''.join(record))
    if not ac_ms['MS_TYPE']:
        ac_ms['MS_TYPE'] = re.findall(r'AC\$ANALYTICAL_CONDITION: MS_TYPE (.+)', ''.join(record))
        Version = 1
    ac_ms['MS_TYPE'] = ac_ms['MS_TYPE'][0].strip() if ac_ms['MS_TYPE'] else None

    if Version == 1:
        ac_ms['IONIZATION'] = re.findall(r'AC\$MASS_SPECTROMETRY: IONIZATION (.+)', ''.join(record))[0].strip()
        ac_ms['ION_MODE'] = re.findall(r'AC\$ANALYTICAL_CONDITION: MODE (.+)', ''.join(record))[0].strip()
    else:
        ac_ms['ION_MODE'] = re.findall(r'AC\$MASS_SPECTROMETRY: ION_MODE (.+)', ''.join(record))[0].strip()

        # Optional AC$MASS_SPECTROMETRY fields
        ac_ms['COLLISION_ENERGY'] = re.findall(r'AC\$MASS_SPECTROMETRY: COLLISION_ENERGY (.+)', ''.join(record))[0].strip()
        # ... (add other optional fields similarly)

        # AC$CHROMATOGRAPHY fields
        ac_lc = {}
        #ac_lc['CAPILLARY_VOLTAGE'] = re.findall(r'AC\$CHROMATOGRAPHY: CAPILLARY_VOLTAGE (.+)', ''.join(record))[0].strip()
        # ... (add other AC$CHROMATOGRAPHY fields similarly)

        # MS$FOCUSED_ION fields
        ms_fi = {}
        ms_fi['BASE_PEAK'] = float(re.findall(r'MS\$FOCUSED_ION: BASE_PEAK (.+)', ''.join(record))[0].strip())
        ms_fi['PRECURSOR_M/Z'] = re.findall(r'MS\$FOCUSED_ION: PRECURSOR_M/Z (.+)', ''.join(record))[0].strip()
        ms_fi['PRECURSOR_TYPE'] = re.findall(r'MS\$FOCUSED_ION: PRECURSOR_TYPE (.+)', ''.join(record))[0].strip()

        if ac_ms['MS_TYPE'] == 'MS2':
            ms_fi['PRECURSOR_M/Z'] = float(ms_fi['PRECURSOR_M/Z'])

    # Handle missing values for ac_ms, ac_lc, ms_fi (simplified)
    for data_dict in [ac_ms, ac_lc, ms_fi]:
        for key in data_dict:
            if not data_dict[key]:
                data_dict[key] = None

    recData['AC$MASS_SPECTROMETRY'] = ac_ms
    recData['AC$CHROMATOGRAPHY'] = ac_lc
    recData['MS$FOCUSED_ION'] = ms_fi

    # Peak extraction
    # PKStart = record.index('PK$PEAK:\n') + 1
    # endslash = record.index('//\n')
    # PKPeak = []
    # for line in record[PKStart:endslash]:
    #     mz, intensity, intrel = map(float, line.strip().split())
    #     PKPeak.append([mz, intensity, intrel])

    # # Annotation extraction (simplified - assuming only one annotation per peak)
    # if readAnnotation:
    #     PKannotationStart = record.index('PK$ANNOTATION:\n') + 1
    #     numpeak = record.index('PK$NUM_PEAK:\n')
    #     if PKannotationStart < numpeak:
    #         for i, line in enumerate(record[PKannotationStart:numpeak]):
    #             mz, formula, formulaCount, mzCalc, dppm = line.strip().split()
    #             PKPeak[i].extend([formula, int(formulaCount), float(mzCalc), float(dppm)])

    # Create DataFrame from PKPeak (using pandas)
    import pandas as pd
    columns = ["mz", "intensity", "intrel"]
    if readAnnotation:
        columns.extend(["formula", "formulaCount", "mzCalc", "dppm"])
    # PKPeak_df = pd.DataFrame(PKPeak, columns=columns)

    # Handle missing values in recData (simplified)
    for key in recData:
        if not recData[key]:
            recData[key] = None

    if recData["CH$SMILES"] == "NA":
        recData["CH$SMILES"] = ""
    if recData["CH$FORMULA"] == "NA":
        recData["CH$FORMULA"] = ""

    # Compose the RmbSpectrum2 object
    sp = RmbSpectrum2()
    # sp = setData(sp, PKPeak_df)  # Set the DataFrame to the data attribute
    sp.info = recData

    # Extract information from sp.info
    sp.info['res'] = sp.info['AC$MASS_SPECTROMETRY'].get('RESOLUTION')
    sp.info['ce'] = sp.info['AC$MASS_SPECTROMETRY'].get('COLLISION_ENERGY')
    sp.info['mode'] = sp.info['AC$MASS_SPECTROMETRY'].get('FRAGMENTATION_MODE')

    return sp
