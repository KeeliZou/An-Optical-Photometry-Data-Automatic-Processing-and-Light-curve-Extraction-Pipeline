from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
import warnings
import glob

warnings.filterwarnings("ignore")

#用于合并BIAS的函数
#raw_bias:原始bias所在路径
#super_bias:超级bias保存路径
def bias_combine(raw_bias,super_bias):
    bias_files = glob.glob(raw_bias)
    bias_files.sort()
    allbias = []

    print("\033[1;31m combining bias ...\033[0m")

    # 合并平场产生超级平场
    for i,ifile in enumerate(bias_files):
        #print("reading bias:",i+1,len(bias_files))
        data = fits.getdata(ifile)
        allbias.append(data)

    print("completed!")
    allbias = np.stack(allbias)
    print(allbias.shape)

    #按行取中位数
    superbias = np.median(allbias,axis=0)
    fits.writeto(super_bias,superbias.astype('float32'),overwrite=True)

#用于产生各个波段的超级平场FLAT
#raw_flat:原始flat所在路径
#super_flat:超级flat保存路径
#super_bias:超级bias保存路径
#band_info:波段信息
def flat_combine(raw_flat,super_flat,super_bias,band_info):
    for j in range(len(band_info)):
        ffiles = glob.glob(raw_flat[j])
        ffiles.sort()
        allflat=[]

        print("\033[1;31m combining "+band_info[j]+" flats...\033[0m")
        superbias = fits.getdata(super_bias)

        for i,ifile in enumerate(ffiles):
            #print("reading "+band_info[j]+" flat:",i + 1,len(ffiles))
            # flat-fielding: subtract bias and then normalize the flat images
            #给FLAT减去BIAS
            data = fits.getdata(ifile) - superbias
            mflat = np.median(data)
            data /= mflat
            #print("median flat:",mflat)
            allflat.append(data)
        allflat = np.stack(allflat)
        print(allflat.shape)
        print("completed!")
        domeflat = np.median(allflat,axis = 0)
        #display the super flat
        #plt.figure(figsize=(8,8))
        #plt.imshow(domeflat,origin='lower')
        #plt.colorbar()
        #plt.title(band_info[j]+" band dome flat drived from dome flats")
        fits.writeto(super_flat[j],domeflat.astype('float32'),overwrite = True)