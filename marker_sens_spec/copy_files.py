out_folder=r'C:\Users\abel\Documents\hydrogenases\marker_sens_spec/filtered'
import os
import shutil
root_dir=r'C:\Users\abel\Documents\hydrogenases\marker_sens_spec'
for fold in os.listdir(root_dir):
    #skip if not a dir

    if os.path.isdir(os.path.join(root_dir,fold)):

        folde=os.path.join(root_dir,fold,"filtered")
        #if fold no exist, skip
        if not os.path.exists(folde):
            continue
        #copy the folder to the out_folder and rename it to the folder name
        shutil.copytree(folde,os.path.join(out_folder,fold))
        print(f"copied {folde} to {os.path.join(out_folder,fold)}")