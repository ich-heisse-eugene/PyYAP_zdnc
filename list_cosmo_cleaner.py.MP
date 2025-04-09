from LA_Cosmic import detCos
import shutil
import os
import logging

import multiprocessing

##################################################################

def cleanout(line, f_out):
    name = line.strip()
    out_name = os.path.splitext(name)[0] + '_CRR.fits'
    detCos(image=name,  out_clean=out_name)
    print(out_name, file=f_out)
    print(out_name)
    logging.info(f"{out_name}")
    print()
    return None

def list_cosmo_cleaner(dir_name, list_name, out_list_name):

    f_out=open(dir_name.joinpath(out_list_name), 'a')

    if multiprocessing.cpu_count() >=2:
        nCPUs = multiprocessing.cpu_count() - 1
    else:
        nCPUs = 1

    with open(dir_name.joinpath(list_name), 'r') as f:

        mp_args = [(curf, f_out) for curf in f]
        with multiprocessing.Pool(processes=nCPUs) as pool:
            pool.starmap_async(cleanout, mp_args, chunksize=nCPUs)
            pool.close()
            pool.join()

    f_out.close()

    print('File names saved in ', dir_name.joinpath(out_list_name))
    return("Cleaned")
