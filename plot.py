import glob
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from astropy.table import Table
from astropy.stats import SigmaClip
from astropy.stats import sigma_clipped_stats
import os


def extract_number(filename):
    # 从文件名中提取数字部分
    number_str = ''.join(filter(str.isdigit, filename))
    return int(number_str)



def auto_plot(band_info,target_star_path,ref_star_cat,check_star_cat,confirmed_star_candidate_path,png_file_path,confirmed_star_path):
    print("\033[1;31m Ploting Images!\033[0m")
    for j in range(len(band_info)):
        
        print("compiling "+band_info[j]+" band images")
        #定义sigmaclip
        sigma_clip = SigmaClip(sigma = 3.0 ,  maxiters = 3)
        
        sfile = target_star_path[j]
        tar_table = Table.read(sfile)

        #参考星星等数组（盛放的是经过一次归一化的数组）
        #此数组是否需要建立掩码看后续情况
        ref_mag = np.empty((len(ref_star_cat[j]),len(tar_table)))
        #参考星星等数组对应mask值
        mask = np.zeros((len(ref_star_cat[j]),len(tar_table)),dtype = bool)
        #合并，将参考星星等数组变成带掩码的数组
        ref_mag = np.ma.masked_array(ref_mag,mask)
        #参考星星等误差均值
        ref_magerr_ave = []
        #目标星星等数组
        tar_mag = np.empty(len(tar_table))
        #目标星星等误差
        tar_magerr = np.empty(len(tar_table))
        #检验星星等误差
        check_magerr = np.empty(len(tar_table))
        
        #mjd数组
        mjd = np.empty(len(tar_table),np.float64)
        #参考星标准差数组
        ref_std_array = np.empty(len(ref_star_cat[j]))
        #虚拟参考星数组
        diff_mag = np.empty(len(tar_table))

        for a in range(len(tar_table)):
            mjd[a] = np.nan
            tar_mag[a] = np.nan
            for b in  range(len(ref_star_cat[j])):
                ref_mag[b,a] = np.nan

        #目标星星等注入（经过初步归一化）
        temp_tar_mag = tar_table['mag']
        tar_mag = temp_tar_mag - np.nanmean(temp_tar_mag)
        tar_magerr = tar_table['magerr']
        mjd = tar_table['mjd']


        #参考星星等注入（经过初步归一化）
        cfile = glob.glob(confirmed_star_candidate_path[j])
        cfile = sorted(cfile,key=extract_number)
        for k in range(len(ref_star_cat[j])):

            ref_table = Table.read(cfile[ref_star_cat[j][k]])
            
            ref_magerr_ave.append(np.nanmean(ref_table['magerr']))

            temp_ref_mag = ref_table['mag']

            ref_mean,ref_median,ref_std = sigma_clipped_stats(temp_ref_mag,sigma=2.0,maxiters=3)

            ref_mag[k] = temp_ref_mag - ref_mean

        #检验星星等注入（经过初步归一化）
        check_table = Table.read(cfile[check_star_cat[j]])
        check_magerr = check_table['magerr']
        temp_check_mag = check_table['mag']
        check_mean,check_median,check_std = sigma_clipped_stats(temp_check_mag,sigma=2.0,maxiters=3)
        check_mag = temp_check_mag - check_mean
        #检验星对应mask值
        mask = np.zeros(len(tar_table),dtype = bool)
        #合并，将检验星星等数组和星等误差数组变成带掩码的数组
        check_mag = np.ma.masked_array(check_mag,mask)
        check_magerr = np.ma.masked_array(check_magerr,mask)


        #对参考星进行sigma-clip操作并计算标准差
        for k in range(len(ref_star_cat[j])):
            temp_data = sigma_clip(ref_mag[k])
            mask_index = np.where(temp_data.mask == True)
            #进行sigma-clip将被clip掉的数据的掩码设置为True
            ref_mag.mask[k][mask_index[0]] = True
            #计算标准差
            ref_std_array[k] = np.nanstd(ref_mag[k])
        
        #对检验星进行sigma-clip操作
        temp_data = sigma_clip(check_mag)
        mask_index = np.where(temp_data.mask == True)
        #进行sigma-clip将被clip掉的数据的掩码设置为True
        check_mag.mask[mask_index[0]] = True
        check_magerr.mask[mask_index[0]] = True
        
        
        #权重函数
        weights = 1.0 / ref_std_array
        #生成虚拟参考星的数据
        temp_diff_mag = np.average(ref_mag,axis=0,weights=weights)
        #归一化虚拟参考星数组
        diff_mag = temp_diff_mag - np.average(temp_diff_mag)

        #进行目标星较差测光操作
        temp_tar_mag_plot = tar_mag - diff_mag
        tar_mag_plot = temp_tar_mag_plot - np.nanmean(temp_tar_mag_plot)
        #进行检验星较差测光操作
        temp_check_mag_plot = check_mag - diff_mag
        check_mag_plot = temp_check_mag_plot - np.nanmean(temp_check_mag_plot)

        #搞点自己喜欢的颜色
        title_color = (34/255,36/255,73/255)
        tar_color = (65/255,76/255,135/255)
        check_color = (175/255,90/255,118/255)
        err_bar_color = (156/255,179/255,212/255)
        
        #绘图模块
        plt.figure(figsize=(16,8),dpi=400)
        #目标星绘图
        plt.errorbar(mjd,tar_mag_plot,yerr = tar_magerr,fmt = '.',color = tar_color,
                     ecolor = err_bar_color,elinewidth=2.5,capsize=4,label = 'Target Star')
        #检验星绘图
        plt.errorbar(mjd,check_mag_plot + 0.5,yerr = check_magerr,fmt = '.',color = check_color,
                     ecolor = err_bar_color,elinewidth=2.5,capsize=4,label = 'Check Star')
        #反转Y轴
        plt.gca().invert_yaxis()
        #添加rms
        temp = np.square(tar_mag_plot)
        temp = np.nanmean(temp)
        target_rms = round(np.sqrt(temp),6)
        temp = np.square(check_mag_plot)
        temp = np.nanmean(temp)
        check_rms = round(np.sqrt(temp),6)

        x_coor = np.median(mjd) - np.min(mjd)
        plt.text(x_coor/2 + np.min(mjd),np.median(tar_mag_plot),'rms =' + str(target_rms),
                 color = title_color,fontsize = 8)
        plt.text(x_coor/2 + np.min(mjd),np.median(check_mag_plot)+0.45,'rms =' + str(check_rms),
                 color = title_color,fontsize = 8)
        #其他杂项
        plt.legend()    
        plt.xlabel('MJD',fontsize = 20)
        plt.ylabel('Delta m')
        plt.title(band_info[j]+' band lightcurve for Target',color = title_color,fontsize=14)
        plt.savefig(png_file_path+band_info[j]+' band lightcurve for Target', bbx_width='tight', dpi=200)
        plt.close()

        tar_err_mean = np.nanmean(tar_magerr)
        tar_std = np.nanstd(tar_mag_plot)
        check_err_mean = np.nanmean(check_magerr)

        check_mean,check_median,check_std = sigma_clipped_stats(check_mag_plot,sigma=3.0,maxiters = 3)

        #存储较差测光信息
        new_header = fits.Header()
        new_header['tar_mag_std'] = tar_std
        new_header['tar_mag_rms'] = target_rms
        new_header['tar_magerr_ave'] = tar_err_mean
        new_header['check_mag_std'] = check_std
        new_header['check_mag_rms'] = check_rms
        new_header['check_magerr_ave'] = check_err_mean

        table = Table()
        table['mjd'] = mjd
        table['tar_mag'] = tar_mag_plot
        table['tar_magerr'] = tar_magerr
        table['check_mag'] = check_mag_plot
        table['check_magerr'] = check_magerr

        bintable_hdu = fits.BinTableHDU(table,new_header)
        hdu1 = fits.HDUList([fits.PrimaryHDU(),bintable_hdu])
        
        catfile = confirmed_star_path+band_info[j]+' band differential photometry result.fits'
        
        hdu1.writeto(catfile,overwrite=True)

        print(" Completed ! ")