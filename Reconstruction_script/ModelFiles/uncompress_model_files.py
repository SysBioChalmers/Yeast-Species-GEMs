# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-


# parse the rhea database
import os    ##for directory
# uncompress
os.system("tar -zxvf xxxx/GitHub/Yeast-Species-GEMs/Reconstruction_script/ModelFiles/*tar.gz")

# compress
os.system("cd GitHub/Yeast-Species-GEMs/Reconstruction_script/ModelFiles/")
os.system("tar -czvf xml.tar.gz xml")
os.system("tar -czvf txt.tar.gz txt")