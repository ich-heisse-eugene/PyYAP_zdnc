Pkg_path = "/Users/sea/Downloads/eshel_zdanice/PyYAP"

import time
from datetime import datetime, date, time
import os
import glob
import shutil
from sys import path, exit
import argparse
import logging
from pathlib import Path

from lister import lister
from trimmer import trimmer
from fixpix import fixpix
from medianer import medianer
from list_subtractor import list_subtractor
from list_dark_subtractor import dark_subtractor
from thar_combiner import thar_combiner
from tracer import order_tracer
from sl_remover import sl_remover
from list_cosmo_cleaner import list_cosmo_cleaner
from extractor import fox
from blz_correct import extract_blz, remove_blz
from thar_reident import thar_auto
from thar_manager import thar_manager
from disp_add import disp_add
from get_sp_resolv import get_sp_resolv
from update_snr import update_snr

## Disable warnings
import warnings
warnings.simplefilter("ignore")
warnings.filterwarnings("ignore")

Pkg_path = Path(Pkg_path)

#####################################################################################
## Parameters
devices = ['eshel_zdnc'] # List of valid devices

#####################################################################################
## Let's get started
def S_EX(conf):
    logging.basicConfig(filename = Path(conf['Path2Data']).joinpath('pyyap_journal.log'), level=logging.INFO, format='%(asctime)s %(message)s')
    logging.getLogger("matplotlib").propagate = False

    start = datetime.now()
    print('Started at:', start.time())

    logging.info("PyYAP - PYthon Yet Another Pipeline\n by Eugene Semenko and Vadim Krushinsky")
    logging.info(f"Started at: {start.time()}")
    logging.info("==========================================")
    logging.info(f"Data in: \'{conf['Path2Data']}\'")
    logging.info(f"Device of origin: {conf['device']}")
    logging.info("Current configuration:")
    print()

    print("==========================================")
    print("PyYAP - PYthon Yet Another Pipeline\n by Eugene Semenko and Vadim Krushinsky")
    print("Reduction of echelle spectra")
    print(f"Data in: \'{conf['Path2Data']}\'")
    print(f"Device of origin: {conf['device']}")
    print("Current configuration:")
    for keys, values in conf.items():
        print(f"{keys} : {values}")
        logging.info(f"{keys} : {values}")

    Path2Data = Path(conf['Path2Data'])
    Path2Raw = Path2Data.joinpath('raw')

    #### Make directory for backing up the raw frames
    if not os.path.exists(Path2Raw):
        os.makedirs(Path2Raw)

    Path2Temp = Path2Data.joinpath('temp')
    #### Make directory for saving of raw frames
    if not os.path.exists(Path2Temp):
        os.makedirs(Path2Temp)

    ## Check directories, create four lists with biases, thars, flats, and objects
    files = lister(Path2Data, Path2Raw, Path2Temp)
    if files != None:
        print ('Lists created')
        logging.info('Lists created')
        print()

   #### Trim overscan and flip spectra
    flip = conf['flip']
    area = list(map(int, conf['area'].strip('][').split(',')))
    print("Trim overscan and flip image")

    ##### Fix CCD cosmetics
    if conf['mask'] != 'None':
        print("Fix cosmetics")
        mask = conf['mask']
        status = fixpix(Path2Data, Path2Temp.joinpath('flat_list.txt'), mask, area, flip)
        print("flat frames ... done")
        logging.info("Cosmetics of flat frames is fixed")
        status = fixpix(Path2Data, Path2Temp.joinpath('thar_list.txt'), mask, area, flip)
        print("ThAr frames fixed")
        logging.info("Cosmetics of ThAr frames is fixed")
        status = fixpix(Path2Data, Path2Temp.joinpath('obj_list.txt'), mask, area, flip)
        print("scientific frames ... done")
        logging.info("Cosmetics of scientific frames is fixed")

    trimmer_data = trimmer(Path2Data, Path2Temp.joinpath('bias_list.txt'), area, flip)
    print("Biases trimmed")
    logging.info("Biases trimmed")

    trimmer_data = trimmer(Path2Data, Path2Temp.joinpath('flat_list.txt'), area, flip)
    print("Flats trimmed")
    logging.info("Flats trimmed")

    trimmer_data = trimmer(Path2Data, Path2Temp.joinpath('thar_list.txt'), area, flip)
    print("ThArs trimmed")
    logging.info("ThArs trimmed")

    trimmer_data = trimmer(Path2Data, Path2Temp.joinpath('obj_list.txt'), area, flip)
    print("Objects trimmed")
    logging.info("Objects trimmed")
    print()

    ## Remove cosmic hits from images
    status = list_cosmo_cleaner(Path2Temp, 'obj_list.txt', 'obj_CRR_cleaned_list.txt')
    print(f"Objects {status}")
    logging.info(f"Objects {status}")
    print()

    ## Create super bias and delete original files
    s_bias_name = conf['s_bias_name']
    print('Create super bias')
    sbias_data = medianer (Path2Data, Path2Temp.joinpath('bias_list.txt'), s_bias_name)
    print(f"Super bias statistic: Mean = {sbias_data[0]:.2f} Median = {sbias_data[1]:.2f} Sigma = {sbias_data[2]:.2f}")
    logging.info(f"Super bias statistic: Mean = {sbias_data[0]:.2f} Median = {sbias_data[1]:.2f} Sigma = {sbias_data[2]:.2f}")
    print()

   #### Subtract super bias from flats
    print('Start flats cleaning')
    status = list_subtractor(Path2Temp.joinpath('flat_list.txt'), Path2Data.joinpath(s_bias_name), 'Bias')
    print (f"Flats: {status}")
    logging.info(f"Flats: {status}")
    print()

   #### Subtract super bias from thars
    print('Start ThArs cleaning')
    status = list_subtractor(Path2Temp.joinpath('thar_list.txt'), Path2Data.joinpath(s_bias_name), 'Bias')
    print(f"ThAr lamp: {status}")
    logging.info(f"ThAr lamp: {status}")
    print()

   #### Subtract super bias from objects
    print('Start objects cleaning')
    status = list_subtractor(Path2Temp.joinpath('obj_CRR_cleaned_list.txt'), Path2Data.joinpath(s_bias_name), 'Bias')
    print(f"Objects: {status}")
    logging.info(f"Objects: {status}")
    print()

    # Work with darks if provided
    if os.path.isfile(Path2Temp.joinpath('dark_list.txt')):
        print('Start subtracting dark current from a series')
        status = dark_subtractor(Path2Data, Path2Temp, ['flat_list.txt', 'thar_list.txt', 'obj_CRR_cleaned_list.txt'], area, flip, s_bias_name)
        print(f"Dark subtraction: {status}")
        logging.info(f"Dark subtraction: {status}")

    #### Create median super flat and delete old flats
    s_flat_name = conf['s_flat_name']
    print('Start to create super flat')
    sflat_data = medianer (Path2Data, Path2Temp.joinpath('flat_list.txt'), s_flat_name)
    print(f"Super flat statistic: Mean = {sflat_data[0]:.2f} Median = {sflat_data[1]:.2f} Sigma = {sflat_data[2]:.2f}")
    print("Super flat created")
    logging.info("Super flat created")
    logging.info(f"Super flat statistic: Mean = {sflat_data[0]:.2f} Median = {sflat_data[1]:.2f} Sigma = {sflat_data[2]:.2f}")
    print()

   ## Create median super ThAr and delete old files
    print('Search nearest ThArs and combine it to super ThArs')
    thars = thar_combiner(Path2Data, Path2Temp.joinpath('thar_list.txt'))
    print(f"Super ThArs saved in {thars}")
    logging.info(f"Super ThArs saved in {thars}")

    slice_half_width = int(conf['slice_half_width'])
    step = int(conf['step'])
    min_height = int(conf['min_height'])
    aperture = float(conf['aperture'])
    view = eval(conf['view'])
    adaptive = eval(conf['adaptive'])

    ### Create list of files for averaging to produce the files with distinctive orders
    if 's_ordim_name' in conf:
        s_ordim_name = conf['s_ordim_name'].rstrip()
    else:
        s_ordim_name = 's_ordim.fits'
    if 's_ordim_method' in conf:      #  Valid methods: 'hybrid', 'flats', 'objects', and 'remap'
        s_ordim_method = conf['s_ordim_method'].rstrip()
    else:
        s_ordim_method = "hybrid"
    print(f"Method: {s_ordim_method}")
    if s_ordim_method != "remap":
        objects_list = []
        if s_ordim_method == "hybrid" or s_ordim_method == "objects":
            with open(Path2Temp.joinpath('obj_CRR_cleaned_list.txt'), 'r') as ff:
                objects_list = ff.read().splitlines()
                ff.close()
            if s_ordim_method == "hybrid":
                objects_list.append(Path2Data.joinpath(s_flat_name))
            with open(Path2Temp.joinpath('ordim_list.txt'), 'w+') as f:
                print(*objects_list, sep='\n', file=f)
                f.close()
            sordim_data = medianer(Path2Data, Path2Temp.joinpath('ordim_list.txt'), Path2Data.joinpath(s_ordim_name))
        if s_ordim_method == 'flats':
            shutil.copy2(Path2Data.joinpath(s_flat_name), Path2Data.joinpath(s_ordim_name))
        if os.path.isfile(Path2Data.joinpath(s_ordim_name)):
            print(f"Master image {s_ordim_name} for the tracer was created using the method '{s_ordim_method}'")
            logging.info(f"Master image {s_ordim_name} for the tracer was created using the method '{s_ordim_method}'")
        else:
            print(f"Error: Failed to create a file {s_ordim_name} for the tracing using the method '{s_ordim_method}'")
            logging.info(f"Error: Failed to create a file {s_ordim_name} for the tracing using the method '{s_ordim_method}'")
            exit(1)

        ## Trace orders
        print("Start tracing orders")
        order_tracer(Path2Data, s_ordim_name, slice_half_width, step, min_height, aperture, adaptive, view)
        shutil.move(os.fspath(Path2Data.joinpath(s_ordim_name)), os.fspath(Path2Temp))
    else:
        print("Using the existing map of orders without modification")
        logging.info("Using the existing map of orders without modification")
        shutil.copy2(Pkg_path.joinpath('devices', conf['device'], 'traces.txt'), Path2Temp)

    # Remove scattered light
    subtract = eval(conf['subtract'])
    if subtract:
        x_order = int(conf['x_order'])
        y_order = int(conf['y_order'])
        ap_file = Path2Temp.joinpath('traces.txt') ##
        sl_remover_data = sl_remover(Path2Data, Path2Temp, s_flat_name, ap_file, step, x_order, y_order, subtract, view)
        print()

        ap_file = Path2Temp.joinpath('traces.txt')
        f_out=open(Path2Temp.joinpath('obj_sl_cleaned_list.txt'), 'a')
        with open(Path2Temp.joinpath('obj_CRR_cleaned_list.txt'), 'r') as f:
            for line in f:
                name = line.strip()
                name = name.split(os.sep)[-1]
                sl_remover_data = sl_remover(Path2Data, Path2Temp, name, ap_file, step, x_order, y_order, subtract, view)
                shutil.move(os.fspath(Path2Data.joinpath(name)), os.fspath(Path2Temp))
                print(sl_remover_data, file=f_out)
                print()
            f.close()
        f_out.close()
        flat_name = Path2Data.joinpath('s_flat_SLR.fits')
        list_name = Path2Temp.joinpath('obj_sl_cleaned_list.txt')
    else:
        flat_name = Path2Data.joinpath('s_flat.fits')
        list_name = Path2Temp.joinpath('obj_CRR_cleaned_list.txt')

    #### Extract objects
    ex_type = conf['ex_type']
    ap_file = Path2Temp.joinpath('traces.txt')
    if ex_type !='FOX':
        status, status_err = fox(flat_name, flat_name, ap_file, ex_type, aperture)
        print(f"Extracted s_flat spectrum saved in {status}")
        logging.info(f"Extracted s_flat spectrum saved in {status}")
    print()

   #### params:
    out_list_name = Path2Temp.joinpath('obj_extracted.txt')
    out_err_list_name = Path2Temp.joinpath('err_extracted.txt')
    f_out=open(out_list_name, 'a')
    fe_out=open(out_err_list_name, 'a')
    with open(list_name, 'r') as f:
        for line in f:
            name = line.strip()
            status, status_err = fox(name, flat_name, ap_file, ex_type, aperture)
            print(f"Extracted spectrum saved in {status}")
            logging.info(f"Extracted spectrum saved in {status}")
            shutil.move(os.fspath(name), os.fspath(Path2Temp))
            print(status, file=f_out)
            print(status_err, file=fe_out)
            print()
    f.close()
    f_out.close()
    fe_out.close()
    os.remove(list_name)
    print()

   #### Extract ThAr
    list_name = Path2Temp.joinpath('s_thar_list.txt')
    out_list_name = Path2Temp.joinpath('s_thar_extracted.txt')
    f_out=open(out_list_name, 'a')
    with open(list_name, 'r') as f:
        for line in f:
            name = line.strip()
            print(f"flat_name = {flat_name}, ap_file = {ap_file}, aperture = {aperture}")
            status, _ = fox(name, flat_name, ap_file, 'APEX', aperture) # Note the fixed method of extraction
            print(f"Extracted spectrum saved in {status}")
            logging.info(f"Extracted spectrum saved in {status}")
            shutil.move(os.fspath(name), os.fspath(Path2Temp))
            print(status, file=f_out)
            print()
    f.close()
    f_out.close()
    os.remove(list_name)
    print()

    # Exit here in case of extraction without wavelength calibration
    if not eval(conf['calibrate']):
        try:
            os.mkdir(Path2Data.joinpath('Reduced'))
        except OSError as e:
            print(f"Exception: {e}")
        finally:
            source = os.listdir(Path2Data)
            for files in source:
                if files.endswith('_ec.fits'):
                    shutil.move(os.fspath(files), os.fspath(Path2Data.joinpath('Reduced')))
                if files.endswith('_err.fits'):
                    shutil.move(os.fspath(files), os.fspath(Path2Data.joinpath('Reduced')))
                if files.endswith('.pdf'):
                    shutil.move(os.fspath(files), os.fspath(Path2Data.joinpath('Reduced')))
                if files.endswith('.log'):
                    shutil.move(os.fspath(files), os.fspath(Path2Data.joinpath('Reduced')))

        if eval(conf['strip']):
            print("Cleaning catalogues...")
            logging.info("Cleaning catalogues...")
            shutil.rmtree(Path2Data.joinpath('temp'),  ignore_errors=True)
            for filepath in glob.iglob(Path2Data.joinpath('*.*')):
                os.remove(filepath)

        end = datetime.now()
        print(f"Finished at: {end.time()}")
        logging.info(f"Finished at: {end.time()}")
        i, d = divmod((end-start).seconds/60, 1)
        print(f"Duration (m:s): {i:3.0f}:{int(d*60)}")
        logging.info(f"Duration (m:s): {i:3.0f}:{int(d*60)}")
        return None
    ### End of the exit point

   #### Search and identify lines in ThAr
   #### params:
    list_name = Path2Temp.joinpath('s_thar_extracted.txt')                  #name of list with extracted thars
    shutil.copy2(Pkg_path.joinpath('devices', conf['device'], conf['device']+'_thar.dat'), Path2Data.joinpath('thar.dat'))                  #copy files to the working directory
    shutil.copy2(Pkg_path.joinpath('devices', conf['device'], conf['device']+'_thar_last.dat'), Path2Data.joinpath('thar_last.dat'))        #copy files to the working directory
    shutil.copy2(Pkg_path.joinpath('devices', conf['device'], conf['device']+'_thar_last.fits'), Path2Data.joinpath('thar_last.fits'))      #copy files to the working directory
    if not os.path.exists(Path2Data.joinpath('Old_Thars')):
        os.makedirs(Path2Data.joinpath('Old_Thars'))

    with open(list_name, 'r') as f:
        anr = float(conf['anr'])
        xord = int(conf['xord'])
        yord = int(conf['yord'])
        for line in f:
            name = line.strip()
            print('Search WL solution for:'+name)
            logging.info("Search WL solution for: {name}")
            thar_auto(Path2Data, name, anr, xord, yord, view)
           ## Save the new WL solution to archive
            dt=datetime.now()
            dt=dt.strftime("%y-%m-%dT%H-%M")
            shutil.copy2(name, Path('Old_Thars/' + dt + '_thar.fits'))
            shutil.copy2(name.replace('_ec', '_err'), Path('Old_Thars/' + dt + '_thar_err.fits'))
            shutil.copy2(os.path.splitext(name)[0] + "_disp.txt", Path('Old_Thars/' + dt + '_thar_disp.txt'))
            shutil.copy2(os.path.splitext(name)[0] + "_features.txt", Path('Old_Thars/' + dt + '_thar_features.txt'))
            print()
    f.close()
    print()

    shutil.move(os.fspath(flat_name), os.fspath(Path2Temp))
    os.remove(Path2Data.joinpath('thar_last.dat'))
    os.remove(Path2Data.joinpath('thar.dat'))
    os.remove(Path2Data.joinpath('thar_last.fits'))
    dt=datetime.now()
    dt=dt.strftime("%y-%m-%dT%H-%M")
    shutil.copy2(Path2Temp.joinpath('traces.txt'), Path('Old_Thars/' + dt + '_traces.txt'))

    ## Search for a nearest calibration for every image and apply the WL solution
    obj_list = Path2Temp.joinpath('obj_extracted.txt')
    obj_err_list = Path2Temp.joinpath('err_extracted.txt')
    thar_list = Path2Temp.joinpath('s_thar_extracted.txt')

    # ThAr
    with open(thar_list, 'r') as f:
        for line in f:
            name = line.strip()
            solution = thar_manager(name, thar_list)
            print(f'{name}: nearest thar is '+str(solution[0].split(os.sep)[-1])+', difference: ' + str(solution[1]) + 'sec')
            logging.info(f'{name}: nearest thar is '+str(solution[0].split(os.sep)[-1])+', difference: ' + str(solution[1]) + 'sec')
            name_cal = disp_add(name, solution[0], view)
            get_sp_resolv(name_cal)
            print()
    f.close()

    # Variances
    with open(obj_err_list, 'r') as f:
        for line in f:
            name = line.strip()
            print(name)
            solution = thar_manager(name, thar_list)
            print(f"Nearest thar is {str(solution[0].split(os.sep)[-1])}, difference: {str(solution[1])} sec")
            logging.info(f"Nearest thar is {str(solution[0].split(os.sep)[-1])}, difference: {str(solution[1])} sec")
            name_cal = disp_add(name, solution[0], view)
            print()
    f.close()

    # Objects
    with open(obj_list, 'r') as f:
        for line in f:
            name = line.strip()
            print(name)
            solution = thar_manager(name, thar_list)
            print(f"{name}: nearest thar is {str(solution[0].split(os.sep)[-1])}, difference: {str(solution[1])} sec")
            logging.info(f"{name}: nearest thar is {str(solution[0].split(os.sep)[-1])}, difference: {str(solution[1])} sec")
            name_cal = disp_add(name, solution[0], view)
            update_snr(name_cal)
            print()
    f.close()

    # Flat and Blaze function
    if ex_type !='FOX':
        name = Path(str(flat_name).replace('.fits', '_ec.fits'))
        solution = thar_manager(name, thar_list)
        print(f"{name}: nearest thar is {str(solution[0].split(os.sep)[-1])}, difference: {str(solution[1])} sec")
        logging.info(f"{name}: nearest thar is {str(solution[0].split(os.sep)[-1])}, difference: {str(solution[1])} sec")
        name_cal = disp_add(name, solution[0], view)
        blz_name = name_cal
        print()

        print("Blaze correction")
        obj_list = Path2Temp.joinpath('obj_extracted.txt')
        with open(obj_list, 'r') as f:
            for line in f:
                name = line.strip()
                in_name = Path(str(name).replace('_ec.fits', '_ec_WCS.fits'))
                out_name = Path(str(name).replace('_ec.fits', '_blz_ec_WCS.fits'))
                status = remove_blz(in_name, blz_name, out_name)
                print(f"Blaze corrected spectrum was saved in {status}")
                logging.info(f"Blaze corrected spectrum was saved in {status}")
        print()

    try:
        if not os.path.exists(Path2Data.joinpath('Reduced')):
            os.mkdir(Path2Data.joinpath('Reduced'))
    except OSError as e:
        print(f"Exception: {e}")
    finally:
        source = os.listdir(Path2Data)
        for files in source:
            if os.path.isfile(files) and not files.endswith('.log'):
                if files.endswith('_WCS.fits'):
                    shutil.copy2(os.fspath(files), os.fspath(Path2Data.joinpath('Reduced')))
                elif files.endswith('.pdf'):
                    shutil.copy2(os.fspath(files), os.fspath(Path2Data.joinpath('Reduced')))
                elif files.endswith('.log'):
                    shutil.copy2(os.fspath(files), os.fspath(Path2Data.joinpath('Reduced')))
                else:
                    shutil.copy2(os.fspath(files), os.fspath(Path2Data.joinpath('temp')))
                os.remove(files)

    # shutil.move(os.fspath(Path2Data.joinpath(s_flat_name)), os.fspath(Path2Temp))
    if eval(conf['strip']):
        print("Cleaning catalogues...")
        logging.info("Cleaning catalogues...")
        shutil.rmtree(Path2Data.joinpath('temp'),  ignore_errors=True)
        shutil.rmtree(Path2Data.joinpath('Old_Thars'),  ignore_errors=True)

    end = datetime.now()
    print(f"Finished at: {end.time()}")
    logging.info(f"Finished at: {end.time()}")
    i, d = divmod((end-start).seconds/60, 1)
    print(f"Duration (m:s): {i:3.0f}:{int(d*60)}")
    logging.info(f"Duration (m:s): {i:3.0f}:{int(d*60)}")
    logging.shutdown()
    shutil.copy2(os.fspath(Path2Data.joinpath('pyyap_journal.log')), os.fspath(Path2Data.joinpath('Reduced')))
    os.remove(os.fspath(Path2Data.joinpath('pyyap_journal.log')))
    return None

## Main programme  ##
if __name__ == "__main__":
    conf = {}
    parser = argparse.ArgumentParser()
    parser.add_argument("cwd", help="Current work directory")
    parser.add_argument("--device", help="Specify the origin of input data", type=str, default="eshel_zdnc")
    parser.add_argument("--extractonly", help="Extraction without calibration", action="store_true")
    parser.add_argument("--sl", help="Turn on/off the subtraction of scattered light [True/False]", default=None)
    parser.add_argument("--method", help="Method of extraction [FOX, APEX, PSFEX]", default=None)
    parser.add_argument("--strip", help="Clean supplement materials [True/False]", default=None)
    parser.add_argument("--view", help="Turn on/of the visualization of the results [True/False]", default=None)
    args = parser.parse_args()

    Data_path = Path(args.cwd.strip())
    conf['Path2Data'] = Path(Data_path)
    conf['device'] = args.device.strip()
    if os.path.isfile(Data_path.joinpath('names.txt')):
        from check_header import fill_headers
        print(f"Correcting FITS headers in {Data_path} \r")
        fill_headers(Data_path.joinpath('names.txt'), conf['device'])
        print()
    if not os.path.isfile(Data_path.joinpath(args.device + '.conf')):
        shutil.copy(Pkg_path.joinpath('devices', args.device, args.device+'.conf'), Data_path)
    with open(Data_path.joinpath(args.device + '.conf')) as f:
        for line in f.readlines():
            if not line.strip().startswith("#"):
                (key, val) = line.split('=')
                conf[key.strip()] = val.strip()

    if args.sl != None and args.sl in ['True', 'False']:
        conf['subtract'] = args.sl
    if args.method != None and args.method in ['FOX', 'APEX', 'PSFEX']:
        conf['ex_type'] = args.method
    if args.strip != None and args.strip in ['True', 'False']:
        conf['strip'] = args.strip
    if args.view != None and args.view in ['True', 'False']:
        conf['view'] = args.view
    if args.extractonly:
        conf['calibrate'] = 'False'
    if conf['mask'].strip() != 'None':
        if not os.path.isfile(Data_path.joinpath(conf['mask'].strip())):
            conf['mask'] = Pkg_path.joinpath('devices', args.device, conf['mask'].strip())

    S_EX(conf)

    exit(0)
