import glob
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
from astropy.time import Time
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.stats import SigmaClip
from astropy.stats import sigma_clipped_stats
from photutils import CircularAperture
from astropy import units as u
import warnings
import os

warnings.filterwarnings("ignore")


#用于判断是否存在TIME-OBS抬头
def has_time_obs_keyword(header):
    if 'TIME-OBS' in header:
        return True
    else:
        return False

#用于进行时间转换
def Time_transfer(file):
    head = fits.getheader(file,ext=1)
    bool_head = has_time_obs_keyword(head)
    datestr = head['DATE-OBS']
    if(bool_head == True):
        timestr = head['TIME-OBS']
        datetime = datestr + 'T' + timestr.strip()
        t = Time(datetime,format='isot',scale='utc')
    else:
        t =Time(datestr,format='isot',scale='utc')
    jd = t.mjd
    return jd

#用于从文件名中提取数字部分
def extract_number(filename):
    number_str = ''.join(filter(str.isdigit, filename))
    return int(number_str)

#用于进行图像对齐、参考星初步选取操作
def star_choose(target_coor,band_info,star_cat_file_path,confirmed_star_path,png_file_path):
    
    #用于进行图像对齐、参考星初步选取操作
    print("\033[1;31m Aligning images and stars!\033[0m")

    #目标坐标转换
    target_coor = SkyCoord(target_coor,unit=(u.hourangle, u.deg))
    coor_ra = target_coor.ra.degree
    coor_dec = target_coor.dec.degree

    #标准图像
    sta_image = []

    for i in range(len(band_info)):
        cfiles = glob.glob(star_cat_file_path[i])
        cfiles.sort()
            
        #计算每一波段中星点数量并做sigma-clip，找出均值
        star_num_array = np.empty(len(cfiles))
        for ii in range(len(cfiles)):
            table = Table.read(cfiles[ii])
            star_num_array[ii] = len(table)

        mean,median,std = sigma_clipped_stats(star_num_array,sigma=3.0)
        #找出最接近中值且大于中值的
        diff_above_median = star_num_array[star_num_array > median] - median
        diff_above_index = np.argmin(diff_above_median)
        diff_above_median = star_num_array[star_num_array > median]
        temp_index = np.where(star_num_array == diff_above_median[diff_above_index])
        star_num_array_index = temp_index[0]


        #读取对应星表作为标准图像
        sta_table = Table.read(cfiles[star_num_array_index[0]])
        sta_head = fits.getheader(cfiles[star_num_array_index[0]],ext=1)
        ref_ra = sta_table['RA']#参考图像中每颗星的RA坐标
        ref_dec = sta_table['DEC']#参考图像中每颗星的DEC坐标
        ni = len(sta_table)#参考图像中的星点数量
        #print(f"ref image star number: {ni}")

        #读取参考图像头部
        rootname,_ = os.path.splitext(cfiles[star_num_array_index[0]])
        ref_head = fits.getheader(rootname[:-4] +'.fits')
        ref_w = WCS(ref_head)

        #将对齐标准图像地址放入标准图像数组
        sta_image.append(rootname[:-4] +'.fits')

        #创建星点矩阵，每一行为一张图像，每一列为同一颗星点
        matrix_mag = np.empty((len(cfiles),ni),np.float32)
        matrix_magerr = np.empty((len(cfiles),ni),np.float32)
        #matrix_time为图像拍摄时间，由于需要精确到秒所以需要float64
        matrix_time = np.empty(len(cfiles),np.float64)
        #matrix_xcenter为每幅图像中星点的x轴坐标
        matrix_xcenter = np.empty((len(cfiles),ni),np.float16)
        #matrix_ycenter为每幅图像中星点的y轴坐标
        matrix_ycenter = np.empty((len(cfiles),ni),np.float16)

        #赤经方向望远镜移动
        ra_movement = np.empty(len(cfiles),np.float64)
        #赤纬方向望远镜移动
        dec_movement = np.empty(len(cfiles),np.float64)

        #赋时+赋空值
        for kk in range(len(cfiles)):
            matrix_time[kk] = Time_transfer(cfiles[kk])
            ra_movement[kk] = np.nan
            dec_movement[kk] = np.nan
            for kkk in range(ni):
                matrix_mag[kk][kkk] = np.nan
                matrix_magerr[kk][kkk] = np.nan

        #将标准图像的星等与星等误差、星点录入指定位置+注入标准图像观测时间
        matrix_mag[star_num_array_index[0]] = sta_table['mag_1']
        matrix_magerr[star_num_array_index[0]] = sta_table['magerr_1']
        matrix_xcenter[star_num_array_index[0]] = sta_table['xcenter']
        matrix_ycenter[star_num_array_index[0]] = sta_table['ycenter']

        #将标准图像中心点赤经赤纬注入
        xcenter = ref_head['NAXIS1']/2
        ycenter = ref_head['NAXIS2']/2
        ra = ref_w.pixel_to_world(xcenter,ycenter).ra.degree
        dec = ref_w.pixel_to_world(xcenter,ycenter).dec.degree
        ra_movement[star_num_array_index[0]] = ra
        dec_movement[star_num_array_index[0]] = dec

        #星点判断阈值
        distance_tolerance = 0.00083
        

        #匹配星点
        for j in range(len(cfiles)):
            #排除参考图像
            if(j == star_num_array_index[0]):
                continue
            else:
                table = Table.read(cfiles[j])
                can_ra = table['RA']
                can_dec = table['DEC']
                can_mag = table['mag_1']
                can_magerr = table['magerr_1']
                can_xcenter = table['xcenter']
                can_ycenter = table['ycenter']

                #读取该星表对应图像，测算图像中心点赤经赤纬
                rootname,_ = os.path.splitext(cfiles[j])
                head = fits.getheader(rootname[:-4]+'.fits')
                w = WCS(head)
                xcenter = head['NAXIS1']/2
                ycenter = head['NAXIS2']/2
                ra = w.pixel_to_world(xcenter,ycenter).ra.degree
                dec = w.pixel_to_world(xcenter,ycenter).dec.degree
                ra_movement[j] = ra
                dec_movement[j] = dec
        
                #寻找最接近目标星点的候选星点
                #计算公式（[delta_RA * cos(DEC)]^2 + [delta_DEC]^2 <= 9 arcsec^2）
                for jj in range(len(ref_ra)):

                    #计算差值
                    delta_ra = np.abs(ref_ra[jj] - can_ra)
                    delta_dec = np.abs(ref_dec[jj] - can_dec)

                    #delta_ra * cos(DEC)
                    ra_diff = delta_ra * np.cos(can_dec * u.degree)

                    #[ra_diff]^2 + [delta_DEC]^2
                    distance = np.square(ra_diff) + np.square(delta_dec)

                    star_indices = np.argmin(distance)

                    #若该星点符合与目标星点赤经赤纬均小于3角秒则认为是同一星点
                    if(distance[star_indices] < np.square(distance_tolerance)):
                        matrix_mag[j][jj] = can_mag[star_indices]
                        matrix_magerr[j][jj] = can_magerr[star_indices]
                        matrix_xcenter[j][jj] = can_xcenter[star_indices]
                        matrix_ycenter[j][jj] = can_ycenter[star_indices]

        #计算星点数量
        valid_data_counts = np.sum(~np.isnan(matrix_mag),axis=0)
        #动态标准，认为其为数据良好的星点
        valid_data_sta = np.average(valid_data_counts)
        valid_data_sta = valid_data_sta + 0.8*(np.max(valid_data_counts)-valid_data_sta)
        valid_data_index = np.where(valid_data_counts > valid_data_sta)

        #进行归一化操作
        #valid_star_mag = np.empty(len(valid_data_index[0]),np.float32)
        #for k in range(len(valid_data_index[0])):
            #temp_mean,temp_median,temp_std = sigma_clipped_stats(matrix_mag[:,valid_data_index[0][k]],sigma=3.0,maxiters = 3)
            #valid_star_mag[k] = temp_mean
        #进行排序得到星点数量较佳且最亮的二十颗星
        #bright_star_index = np.argsort(valid_star_mag)
        #for j in range(len(cfiles)):
            #final_mean,final_median,final_std = sigma_clipped_stats(matrix_mag[j,valid_data_index[bright_star_index[:-20]]],sigma = 3.0,maxiters = 3)


        
        
        #在计算标准差之前先对星点进行归一化操作
        
        #从数据点良好的星点中挑选标准差小的星点
        #此操作用来防止参考星是变星
        valid_star_std = np.empty(len(valid_data_index[0]),np.float32)
        for k in range(len(valid_data_index[0])):
            valid_star_std[k] = np.nanstd(matrix_mag[:,valid_data_index[0][k]])
        std_lim = np.average(valid_star_std)
        valid_std_index = np.where(valid_star_std < std_lim)
        #数据点多且标准差小的参考星星点索引
        good_star_index = valid_data_index[0][valid_std_index[0]]

        
        #参考星近邻星点阈值
        ref_tolerance = 0.00139

        ra_temp = sta_table['RA']
        dec_temp = sta_table['DEC']

        #存储潜在参考星
        for a in range(len(good_star_index)):
            
            ra = sta_table[good_star_index[a]]['RA']
            dec = sta_table[good_star_index[a]]['DEC']

            #用于剔除临近有其他星点的参考星
        
            ra_dis = np.abs(ra - ra_temp)
            ra_dis = ra_dis * np.cos(dec * u.degree)
            dec_dis = np.abs(dec - dec_temp)
            ra_search_index = np.where(ra_dis < ref_tolerance)
            dec_search_index = np.where(dec_dis < ref_tolerance)
            ref_index = np.intersect1d(ra_search_index,dec_search_index)

            #说明只有潜在参考星本身
            if(len(ref_index) == 1):

                #判断是否为目标星（用于剔除目标星是系外行星宿主恒星这种极端情况）
                temp_ra_dis = np.abs(ra - coor_ra)
                temp_ra_dis = temp_ra_dis * np.cos(coor_dec * u.degree)
                temp_dec_dis = np.abs(dec - coor_dec)
                if(temp_ra_dis > ref_tolerance and temp_dec_dis > ref_tolerance):
                    new_head = fits.Header()
                    new_head['SIMPLE'] = True
                    new_head['BITPIX'] = -32
                    new_head['RA'] = ra
                    new_head['DEC'] = dec

                    #sky = SkyCoord(ra = ra * u.deg,dec = dec * u.deg,frame='icrs')
                    #xcenter,ycenter = ref_w.world_to_pixel(sky)

                    new_head['XCENTER'] = matrix_xcenter[star_num_array_index[0],good_star_index[a]]
                    new_head['YCENTER'] = matrix_ycenter[star_num_array_index[0],good_star_index[a]]

                    #写入数据
                    new_table = Table()
                    new_table['mjd'] = matrix_time
                    new_table['mag'] = matrix_mag[:,good_star_index[a]]
                    new_table['magerr'] = matrix_magerr[:,good_star_index[a]]
                    new_table['xcenter'] = matrix_xcenter[:,good_star_index[a]]
                    new_table['ycenter'] = matrix_ycenter[:,good_star_index[a]]

                    bintable_hdu = fits.BinTableHDU(new_table,new_head)
                    hdu1 = fits.HDUList([fits.PrimaryHDU(),bintable_hdu])
                    #存储参考星
                    hdu1.writeto(confirmed_star_path+band_info[i]+'_band_confirmed_star_'+ str(a+1)+'.fits',overwrite=True)
                else:
                    #说明此源为目标星
                    continue
            else:
                #说明此星点周围5角秒内存在其他星点，抛弃此参考星
                continue
        
        
        
        #用于从table中读取目标星的星等、星等误差保存
        temp_delta_ra = np.abs(coor_ra - sta_table['RA'])

        temp_delta_ra = temp_delta_ra * np.cos(coor_dec * u.degree)

        temp_delta_dec = np.abs(coor_dec - sta_table['DEC'])

        temp_distance = np.square(temp_delta_ra) + np.square(temp_delta_dec)

        distance_indices = np.argmin(temp_distance)

        #此处存在一个安全隐患 即：倘若5角秒内有多个目标会报错
        #此算法被抛弃
        #tar_ra_index = np.where(temp_ra < tar_ra_tolerance)
        #tar_dec_index = np.where(temp_dec < tar_dec_tolerance)
        #tar_index = np.int(np.intersect1d(tar_ra_index,tar_dec_index))

        if(temp_distance[distance_indices] < np.square(distance_tolerance)):


            tar_ra = sta_table[distance_indices]['RA']
            tar_dec = sta_table[distance_indices]['DEC']

            new_head = fits.Header()
            new_head['SIMPLE'] = True
            new_head['BITPIX'] = -32
            new_head['RA'] = tar_ra
            new_head['DEC'] = tar_dec

            xcenter = sta_table[distance_indices]['xcenter']
            ycenter = sta_table[distance_indices]['ycenter']

            new_head['XCENTER'] = np.float32(xcenter)
            new_head['YCENTER'] = np.float32(ycenter)

            #写入数据
            new_table = Table()
            new_table['mjd'] = matrix_time
            new_table['mag'] = matrix_mag[:,distance_indices]
            new_table['magerr'] = matrix_magerr[:,distance_indices]

            bintable_hdu = fits.BinTableHDU(new_table,new_head)
            hdu1 = fits.HDUList([fits.PrimaryHDU(),bintable_hdu])
            #存储目标星
            hdu1.writeto(confirmed_star_path+band_info[i]+'_band_target_star'+'.fits',overwrite=True)
        
        else:
            print('出错！该坐标附近未找到可识别的参考星！')

        # 图像位移
        plt.figure(figsize=(8,6),dpi=100)
        plt.plot(matrix_time,ra_movement)
        plt.title('Target '+band_info[i]+' band RA movement',fontsize=20)
        plt.xlabel('MJD',fontsize = 15)
        plt.ylabel('RA MOVEMENT',fontsize = 15)
        #plt.savefig(png_file_path+band_info[i]+' band RA movement.png', bbx_width='tight',dpi=50)
        #论文使用这个
        plt.savefig(png_file_path+band_info[i]+' band RA movement.eps', bbx_width='tight',dpi=100)
        plt.close()

        plt.figure(figsize=(8,6),dpi=100)
        plt.plot(matrix_time,dec_movement)
        plt.title('Target '+band_info[i]+' band DEC movement',fontsize=20)
        plt.xlabel('MJD',fontsize = 15)
        plt.ylabel('DEC MOVEMENT',fontsize = 15)
        #plt.savefig(png_file_path+band_info[i]+' band DEC movement.png', bbx_width='tight', dpi=50)
        #论文使用这个
        plt.savefig(png_file_path+band_info[i]+' band DEC movement.eps', bbx_width='tight', dpi=100)
        plt.close()

        print(" Completed ! ")

    return sta_image
        



#用于进行参考星进一步选取操作
def refer_chooes(band_info,confirmed_star_candidate_path,target_star_path,sta_image,png_file_path):
    
    print("\033[1;31m Selecting candidate stars!\033[0m")

    #参考星索引列表
    ref_star_cat = []
    #检验星索引列表
    check_star_cat = []
    for j in range(len(band_info)):

    ######计算区域######
        
        #读取目标星HEADER信息
        target_head = fits.getheader(target_star_path[j],ext=1)
        RA = target_head['RA']
        DEC = target_head['DEC']

        #读取目标星table信息
        target_table = Table.read(target_star_path[j])
        tar_mag = target_table['mag']
        tar_mag_ave = np.nanmean(tar_mag)


        #读取参考星列表
        cfile = glob.glob(confirmed_star_candidate_path[j])
        cfile = sorted(cfile,key=extract_number)

        #以(x,y,RA,DEC)格式存储候选星
        candidates = np.empty((len(cfile),4))
        
        #各参考星各图像星等标准差
        #ref_magerr = np.empty((len(cfile),len(tar_mag)))
        #这个就不用了（一般来说星等越亮测光误差也越小）

        #各参考星各图像星等（最后一列为目标星）
        ref_mag = np.empty((len(cfile)+1,len(tar_mag)))
        #各参考星校正后的星等平均值
        ref_mag_ave = np.empty(len(cfile))
        #为星等均值数组设置掩码
        mask = np.zeros(len(ref_mag_ave),dtype = bool)
        ref_mag_ave = np.ma.masked_array(ref_mag_ave,mask)
        #各参考星校正后星等标准差
        ref_mag_std  = np.empty(len(cfile))
        #为星等均值数组设置掩码
        mask = np.zeros(len(ref_mag_std),dtype = bool)
        ref_mag_std = np.ma.masked_array(ref_mag_std,mask)


        #读取参考星FITS文件信息
        for i in range(len(cfile)):
            #读取头部
            head = fits.getheader(cfile[i],ext=1)
            #以(x,y,RA,DEC)格式存储候选星
            candidates[i,0] = head['XCENTER']
            candidates[i,1] = head['YCENTER']
            candidates[i,2] = head['RA']
            candidates[i,3] = head['DEC']

            #读取数据
            table = Table.read(cfile[i])
        
            #将参考星星等赋值给ref_mag
            ref_mag[i] = table['mag']
        
        #最后一列数据赋值
        ref_mag[len(cfile)] = tar_mag
        


        ####此函数用于找出时间最早的拥有所有星点的图像序号
        # 使用 np.isnan 函数检查矩阵中的 NaN 值
        nan_mask = np.isnan(ref_mag)

        # 使用 np.any 函数检查哪些列不存在 NaN 值
        non_nan_columns = ~np.any(nan_mask, axis=0)

        # 获取不存在 NaN 值的列索引
        column_indices = np.where(non_nan_columns)[0]
        
        #此数组用于盛放参考星全的星点均值
        temp_mag_ave = np.empty(len(column_indices))


        #计算其中最亮的图像，将其作为基准图像
        #为1直接不进行循环
        if(len(column_indices) == 1):
            bright_image_index = column_indices[0]
            continue
        else:
            for k in range(len(column_indices)):
                temp_mean,temp_median,temp_std = sigma_clipped_stats(ref_mag[:,column_indices[k]],sigma=3.0)
                temp_mag_ave[k] = temp_mean
            temp_index = np.argmin(temp_mag_ave)
            bright_image_index = column_indices[temp_index]

        #检索最亮图像中的最亮的二十颗星
        bright_star_index = np.argsort(ref_mag[:,bright_image_index])

        #将每一帧数据归一化到第一幅图像整体水平
        for k in range(len(tar_mag)):
            if(k == bright_image_index):
                continue
            else:
                temp_ref_mag = ref_mag[bright_star_index[:-20],bright_image_index] - ref_mag[bright_star_index[:-20],k]
                ref_mean,ref_median,ref_std = sigma_clipped_stats(temp_ref_mag,sigma=3.0)
                ref_mag[:,k] = ref_mag[:,k] + ref_mean

        
        #将定标后的星等重新写入
        for i in range(len(cfile)):
            head = fits.getheader(cfile[i],ext=1)
            table = Table.read(cfile[i])
            #只改变星等一列，其他都不改变
            table['mag'] = ref_mag[i]
            new_header = fits.Header()
            bintable_hdu = fits.BinTableHDU(table,new_header)
            hdu1 = fits.HDUList([fits.PrimaryHDU(),bintable_hdu])
            hdu1.writeto(cfile[i],overwrite=True)
        
        #修改target_star
        head = fits.getheader(target_star_path[j],ext=1)
        table = Table.read(target_star_path[j])
        #只改变星等一列
        table['mag'] = ref_mag[len(cfile)]
        new_header = fits.Header()
        bintable_hdu = fits.BinTableHDU(table,new_header)
        hdu1 = fits.HDUList([fits.PrimaryHDU(),bintable_hdu])
        hdu1.writeto(target_star_path[j],overwrite=True)
        
        
        
    #####绘图区域#########

        #记得传入这个路径
        data = fits.getdata(sta_image[j])

        #潜在参考星
        positions = np.transpose((candidates[:,0],candidates[:,1]))
        apertures = CircularAperture(positions, r = 15.)
        plt.figure(figsize = (8,8))
        #plt.figure(figsize = (16,16))
        plt.imshow(data,cmap = 'Greys_r',origin = 'lower',vmin = np.median(data)-0.5*np.std(data),
                    vmax = np.median(data)+1.8*np.std(data),interpolation = 'nearest')
        apertures.plot(color = 'red',lw = 1.5,alpha = 0.5)
        plt.title(band_info[j]+' band reference stars candidates',fontsize=20)
        for ii in range (len(positions)):
            plt.text(positions[ii][0],positions[ii][1],
                    str((f'{candidates[ii,2]:.5f}',f'{candidates[ii,3]:.5f}')),color='red')
        


        #目标星
        positions = np.transpose((target_head['XCENTER'],target_head['YCENTER']))
        apertures = CircularAperture(positions, r = 15.)
        apertures.plot(color = 'yellow',lw = 1.5,alpha = 0.5)
        plt.text(target_head['XCENTER'],target_head['YCENTER'],str((f'{RA:.5f}',f'{DEC:.5f}')),color='yellow')

        #存储图像
        #plt.savefig(png_file_path+band_info[j]+' band reference stars candidates.png', bbx_width='tight', dpi=600)
        plt.savefig(png_file_path+band_info[j]+' band reference stars candidates.eps', bbx_width='tight', dpi=100)
        plt.close()


    ######计算区域######

        #计算潜在参考星星等均值与误差
        for i in range(len(cfile)-1):
            ref_mean,ref_median,ref_std = sigma_clipped_stats(ref_mag[i],sigma=3.0)
            ref_mag_ave[i] = ref_mean
            ref_mag_std[i] = ref_std
        
        
        ####挑选检验星####
        temp_ref_mag = np.abs(tar_mag_ave - ref_mag_ave)
        #找到所有与目标星星等均值小于一等的参考星（可以暗于也可以亮于）
        candidates_index = np.where(temp_ref_mag <= 1)
        #挑选其中标准差最小者
        temp_check_star_index = np.argmin(ref_mag_std[candidates_index[0]])
        
        check_star_index = candidates_index[0][temp_check_star_index]
        
        
        #将两个数组中检验星对应的掩码设置为True,不参与后续的参考星选取
        ref_mag_std.mask[check_star_index] = True
        ref_mag_ave.mask[check_star_index] = True

        #####选取参考星#####

        #取校正后最亮的20颗星
        mag_sorted_index = np.argsort(ref_mag_ave)
        #在这二十颗星中取标准差最小的10颗星
        std_sorted_index = np.argsort(ref_mag_std[mag_sorted_index[:20]])
        #将这些星作为检验星
        ref_index = std_sorted_index[:10]

    #####绘图区域#########

        #绘制自动选择后的图
        #plt.figure(figsize = (16,16))
        plt.figure(figsize = (8,8))
        plt.imshow(data,cmap = 'Greys_r',origin = 'lower',vmin = np.median(data)-0.5*np.std(data),
                    vmax = np.median(data)+1.8*np.std(data),interpolation = 'nearest')
        plt.title(band_info[j]+' band reference stars(AutoChoose)',fontsize=20)

        #目标星绘图
        positions = np.transpose((target_head['XCENTER'],target_head['YCENTER']))
        apertures = CircularAperture(positions, r = 15.)
        apertures.plot(color = 'yellow',lw = 1.5,alpha = 0.5)
        plt.text(target_head['XCENTER'],target_head['YCENTER'],str((f'{RA:.5f}',f'{DEC:.5f}')),color='yellow')

        # 检验星绘图
        positions = np.transpose((candidates[check_star_index,0],candidates[check_star_index,1]))
        apertures = CircularAperture(positions, r = 15.)
        apertures.plot(color = 'blue',lw = 1.5,alpha = 0.5)
        plt.text(candidates[check_star_index,0],candidates[check_star_index,1],str((f'{RA:.5f}',f'{DEC:.5f}')),color='blue')
        
        #自动选取的10颗参考星
        for ii in range (len(ref_index)):
            positions = np.transpose((candidates[ref_index[ii],0],candidates[ref_index[ii],1]))
            apertures = CircularAperture(positions, r = 15.)
            apertures.plot(color = 'red',lw = 1.5,alpha = 0.5)
            plt.text(candidates[ref_index[ii]][0],candidates[ref_index[ii]][1],
                        f"No: {ii+1}, Index:{ref_index[ii]}",color='red')
        
        #存储图像
        #plt.savefig(png_file_path+band_info[j]+' band reference stars(AutoChoose).png', bbx_width='tight', dpi=600)
        plt.savefig(png_file_path+band_info[j]+' band reference stars(AutoChoose).eps', bbx_width='tight', dpi=100)
        plt.close()


        #添加参考星索引
        ref_star_cat.append(ref_index)
        #添加检验星索引
        check_star_cat.append(check_star_index)

        print(" Completed ! ")

    return ref_star_cat,check_star_cat