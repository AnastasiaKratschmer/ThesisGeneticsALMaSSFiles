#run first:
#pip install pandas
#pip install os.path2
#pip install subprocess.run
#pip install pytest-shutil



import pandas as pd
import numpy as np
import os, os.path
import subprocess
import shutil

#loll

#from common.readcsv import readcsv
#from common.gis2almass import make_farms_df
#from common.gis2almass import make_polyref_df
#from common.ascii2lsb import ascii2lsb
#from common.utm2LatLng import utmToLatLng

print("I am running")
landscapes=["AKK24_LS01_amNS","AKK24_LS02_amNS","AKK24_LS03_amEW","AKK24_LS04_amEW","AKK24_LS05_amEW"]  #The "real" landscapes!
#landscapes=["AKK24_LS06_nm","AKK24_LS07_nm","AKK24_LS08_nm","AKK24_LS09_nm","AKK24_LS10_nm"] #nm landscapes
#landscapes=["AKK24_LS06_fmNS","AKK24_LS07_fmNS","AKK24_LS08_fmNS","AKK24_LS09_fmNS","AKK24_LS10_fmNS"] #fake NS landscapes
#landscapes=["AKK24_LS06_fmEW","AKK24_LS07_fmEW","AKK24_LS08_fmEW","AKK24_LS09_fmEW","AKK24_LS10_fmEW"] #fake EW landscapes
#landscapes=["AKK24_LS06_frNS","AKK24_LS07_frNS","AKK24_LS08_frNS","AKK24_LS09_frNS","AKK24_LS10_frNS"] #"fake real" NS landscapes
#landscapes=["AKK24_LS06_fmEW","AKK24_LS07_fmEW","AKK24_LS08_fmEW","AKK24_LS09_fmEW","AKK24_LS10_fmEW"] #"fake real" EW landscapes
#landscapes=["AKK24_LS01_amNS"]
weather_dict = {
    "01": "LS01_era5.pre",
    "02": "LS02_era5.pre",
    "03": "LS03_era5.pre",
    "04": "LS04_era5.pre",
    "05":"LS05_era5.pre",
    "06":"LS06_era5.pre",
    "07":"LS07_era5.pre",
    "08":"LS08_era5.pre",
    "09":"LS09_era5.pre",
    "10":"LS10_era5.pre"
}

running_folder = "./vole_debugAM/"
pathout= "Countries/test_DK/"

for i in range(0,10):
    for landscape_name in landscapes:
        
        print(landscape_name, " ... started", "iteration", str(i))

        landscape_id=landscape_name[8:10]
        weather_file_string=weather_dict[landscape_id]


        polyref_filename_short =  landscape_name + "_LEbCh25W4_polyref" + ".txt"
        polyref_filename = running_folder+pathout + polyref_filename_short

        farmref_filename_short =  landscape_name + "_LEbCh25W4_farmref" + ".txt"
        farmref_filename = running_folder+pathout + farmref_filename_short

        raster_lsb_filename_short = landscape_name +"_LEbCh25W4"+ ".lsb"
        raster_lsb_filename = running_folder+pathout + raster_lsb_filename_short


        lsb_ALMaSSrun_src = os.path.join(pathout, raster_lsb_filename_short)
        pref_ALMaSSrun_src = os.path.join(pathout, polyref_filename_short)
        fref_ALMaSSrun_src = os.path.join(pathout, farmref_filename_short)
        weather_ALMaSSrun_scr=os.path.join(pathout,weather_file_string)

        lsb_exists = os.path.exists(lsb_ALMaSSrun_src)
        pref_exists = os.path.exists(pref_ALMaSSrun_src)
        fref_exists = os.path.exists(fref_ALMaSSrun_src)

        print("lsb_ALMaSSrun_src exists: ", lsb_ALMaSSrun_src,":",lsb_exists)
        print("pref_ALMaSSrun_src exists:", pref_exists)
        print("fref_ALMaSSrun_src exists:", fref_exists)

        with open(os.path.join(running_folder, 'TIALMaSSConfig.cfg'), 'r+') as file:
            print("I openede the config file!")
            lines = file.readlines()
            file.seek(0)
            file.truncate()
            file.writelines(lines[:-5])

        cfg_file = open(os.path.join(running_folder, 'TIALMaSSConfig.cfg'), 'a')
        cfg_file.write('MAP_POLY_FILE (string) = ' +'\"'+pref_ALMaSSrun_src+'\"'+'\n')
        cfg_file.write('MAP_FARMREF_FILE (string) = ' +'\"'+fref_ALMaSSrun_src+'\"'+'\n')
        cfg_file.write('MAP_MAP_FILE (string) = ' +'\"'+lsb_ALMaSSrun_src+'\"'+'\n')
        cfg_file.write('LANDSCAPE_INFO (string) = ' +'\"'+lsb_ALMaSSrun_src+'\"'+'\n')
        cfg_file.write('MAP_WEATHER_FILE (string) = ' +'\"'+weather_ALMaSSrun_scr+'\"'+'\n')
        print("I wrote to the config file ")
        
        cfg_file.close()
        
        current_pro=subprocess.Popen(['./almass_cmd'], cwd=running_folder)

        current_pro.wait()
        print("\n \n \n job stopped/ended \n \n \n ")

        pathoutput="./output/"+landscape_name+"/rep_"+str(i)

        if not os.path.exists(pathoutput):
            os.makedirs(pathoutput)

        source="./vole_debugAM/DistanceAndGeneticDifference.txt"
        shutil.move(source, pathoutput)

        source="./vole_debugAM/FIS_output_dataset_style.txt"
        shutil.move(source, pathoutput)

        source="./vole_debugAM/FST_output_matrix_style.txt"
        shutil.move(source, pathoutput)

        source="./vole_debugAM/FST_output_dataset_style.txt"
        shutil.move(source,pathoutput)

        source="./vole_debugAM/Heterozyg_output.txt"
        shutil.move(source,pathoutput)

        source="./vole_debugAM/Probe.res"
        shutil.move(source, pathoutput)

        source="./vole_debugAM/PopulationInQuadrants.txt"
        shutil.move(source, pathoutput)

        source="./vole_debugAM/GeneticsAndPosition1.txt"
        shutil.move(source, pathoutput)

        source="./vole_debugAM/GeneticsAndPosition6.txt"
        shutil.move(source, pathoutput)

        source="./vole_debugAM/GeneticsAndPosition11.txt"
        shutil.move(source, pathoutput)

        source="./vole_debugAM/GeneticsAndPosition16.txt"
        shutil.move(source, pathoutput)

        source="./vole_debugAM/GeneticsAndPosition21.txt"
        shutil.move(source, pathoutput)

        source="./vole_debugAM/GeneticsAndPosition26.txt"
        shutil.move(source, pathoutput)

        source="./vole_debugAM/GeneticsAndPosition31.txt"
        shutil.move(source, pathoutput)
