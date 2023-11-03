from astropy.table import Table
from astropy.io import fits
from astroquery.astrometry_net import AstrometryNet
from astropy import units as u
from astropy.coordinates import SkyCoord
import photutils as pht
import warnings
from astropy.stats import SigmaClip
from astropy.stats import sigma_clipped_stats
from photutils.utils import calc_total_error
from astropy.wcs import WCS
import numpy as np
import glob
import os
warnings.filterwarnings("ignore")

#此操作同时进行背景天光估计与找星工作
#step: background_estimate->star_find->mark_star_position



    
#用于判断是否FITS头部存在TIME-OBS关键词
def has_time_obs_keyword(header):
    if 'TIME-OBS' in header:
        return True
    else:
        return False


#估计天光背景，评估哪些图像质量可能不佳不参加天测操作
def bkg_evaluate(band_info,cut_data_file_path):
    #向天测函数传递的cfiles
    print("\033[1;31m Evaluating images background!\033[0m")
    cfiles = []
    for j in range(len(band_info)):
        
        print("compiling "+band_info[j]+" band images")
        sfiles = glob.glob(cut_data_file_path[j])
        sfiles.sort()

        #申请用于存贮background均值的数组
        bkg_median = np.empty(len(sfiles))
        mask = np.zeros(len(sfiles),dtype = bool)
        bkg_median = np.ma.masked_array(bkg_median,mask)

        sigmaclip = SigmaClip(sigma = 2.5,maxiters = 3)

        for i,ifile in enumerate(sfiles):
            data = fits.getdata(ifile)
            mask = pht.make_source_mask(data,nsigma = 5, npixels = 10, dilate_size = 11)
            bkg_estimator = pht.SExtractorBackground()
            bkg = pht.Background2D(data,(64,64),mask = mask,filter_size = (3,3),sigma_clip = SigmaClip(sigma = 3.),
                                    bkg_estimator = bkg_estimator)
            bkg_median[i] = bkg.background_median

            progress = (i+1) / len(sfiles) * 100
            print(f"working process: [{int(progress):3d}%],",end='\r')

        temp_bkg_median = sigmaclip(bkg_median)
        mask_index = np.where(temp_bkg_median.mask == True)
        bkg_median.mask[mask_index[0]] = True
        cfiles_index = np.where(bkg_median.mask == False)
        
        temp_cfiles = []
        for ii in range(len(cfiles_index[0])):
            temp_cfiles.append(sfiles[cfiles_index[0][ii]])
        cfiles.append(temp_cfiles)
        print(" Completed ! ")
    return(cfiles)


#上传Astrometry.net完成解算，形成天球坐标
def ast_solve(band_info,cfiles,fwhm_initial,ast_sigma):
    error_image1 = []
    print("\033[1;31m Requesting Astrometry.net!\033[0m")
    for ii in range(len(band_info)):
        print("compiling "+band_info[ii]+" band images")
        
        #直接使用传入的cfiles,经过背景天光估计的
        for i,ifile in enumerate(cfiles[ii]):
            rootname,_ = os.path.splitext(ifile)
            catfile = rootname + '-ast.fits'
            data = fits.getdata(ifile)
            head = fits.getheader(ifile)
            #用于判断是否存在TIME-OBS关键词
            #对head做坐标转换
            coordinate = head['RA']+head['DEC']
            c = SkyCoord(coordinate,unit=(u.hourangle,u.deg))
            #先进行写入，更改head信息
            ra = c.ra.degree
            dec = c.dec.degree
            head['RA'] = ra
            head['DEC'] = dec
            fits.writeto(ifile,data,header=head,overwrite=True)
            print("Finding stars in"+rootname)
            #进行天测操作
            for attempt in range(2):
                try:
                    ast = AstrometryNet()
                    ast.api_key = 'ujgxlgcypxbylyik'
                    #30倍sigma，只找亮星
                    print(band_info[ii]+'波段第'+str(i+1)+'/'+str(len(cfiles[ii]))+'张')
                    wcs_header = ast.solve_from_image(ifile,ra_dec_units = 'degree',fwhm = fwhm_initial,detect_threshold = ast_sigma, solve_timeout=40)
                    if(wcs_header == {}):
                        print('解算异常！头文件为空，抛弃此图像！' + rootname)
                        error_image1.append(rootname)
                        #此处也需要直接抛出异常
                        break
                    else:
                        #写入新地址
                        fits.writeto(catfile,data,header=wcs_header,overwrite=True)
                        #一次解算成功即退出循环
                        break                    
                
                #解算超时或无法解算处理
                except Exception as e:
                    #判断是否是第一次解算失败
                    if(attempt == 0):
                        print('\033[1;31m 解算超时！重新提交进行解算 \033[0m')
                        #第一次解算失败，再次提交
                        continue
                    else:
                        print('\033[1;31m 解算超时！此图像无法解算！ \033[0m')
                        error_image1.append(rootname)
    
    print(" Completed ! ")

    return error_image1

#用于计算增益与读出噪声
def gain_caculate(band_info,biasfile1,biasfile2,flatfile1,flatfile2):
    print("\033[1;31m Calculating Gain and ReadoutNoise!\033[0m")
    gain = []
    rdnoise = []
    for i in range(len(band_info)):
        bias1 = fits.getdata(biasfile1)
        bias2 = fits.getdata(biasfile2)
        flat1 = fits.getdata(flatfile1[i])
        flat2 = fits.getdata(flatfile2[i])

        mean_flat1 = np.median(flat1)
        mean_flat2 = np.median(flat2)
        mean_bias1 = np.median(bias1)
        mean_bias2 = np.median(bias2)

        _,_,std_biasdiff = sigma_clipped_stats(bias1 - bias2,sigma = 2.0,maxiters = 2)
        _,_,std_flatdiff = sigma_clipped_stats(flat1 - flat2,sigma = 3.0,maxiters = 2)
        #增益计算
        gain.append(((mean_flat1 + mean_flat2) - (mean_bias1 + mean_bias2))/((std_flatdiff**2 - std_biasdiff**2)))
        rdnoise.append(gain[i] * std_biasdiff / np.sqrt(2))
    return gain,rdnoise


#找星函数
def star_find(band_info,ast_data_file_path,gain,fwhm_initial,star_num):
    print("\033[1;31m Finding stars in images!\033[0m")

    error_image2 = []
    
    for ii in range(len(band_info)):
        print("compiling "+band_info[ii]+" band images")
        cfiles = glob.glob(ast_data_file_path[ii])
        cfiles.sort()


        #此操作可以保证只对天测成功的图像进行找星
        for i,ifile in enumerate(cfiles):
            rootname,_ = os.path.splitext(ifile)
            rootname = rootname[:-4] + '.fits'
            data = fits.getdata(ifile)
            raw_head = fits.getheader(rootname)
            ast_head = fits.getheader(ifile)
            #创建WCS
            w = WCS(ast_head)
            catfile = ifile[:-5] +'-cat.fits'

            #判断此FITS头是否含有TIME-OBS关键词
            bool_head = has_time_obs_keyword(raw_head)
            #print("Finding stars in"+rootname)

            #测量图像背景及测光误差计算
            mask = pht.make_source_mask(data,nsigma = 5, npixels = 10, dilate_size = 11)
            bkg_estimator = pht.SExtractorBackground()
            bkg = pht.Background2D(data,(64,64),mask = mask,filter_size = (3,3),sigma_clip = SigmaClip(sigma = 3.),
                                    bkg_estimator = bkg_estimator)
            error = calc_total_error(data - bkg.background,bkg.background_rms, effective_gain= gain[ii])

            #使用固定的半高全宽进行找星、以判断真正的半高全宽值
            mean,median,std = sigma_clipped_stats(data,sigma=3.0)
            daofind = pht.IRAFStarFinder(fwhm = fwhm_initial,minsep_fwhm = 3.0, threshold = 8*std,exclude_border = True,
                                         sharplo = 0.5,sharphi = 2.0, roundlo = 0.0,roundhi = 0.7)
            sources = daofind(data - bkg.background)

            #提取初步半高全宽
            fwhm_raw = sources['fwhm']
            fwhm_mean,fwhm_median,fwhm_std = sigma_clipped_stats(fwhm_raw,sigma=3.0)

            #再次找星
            daofind = pht.IRAFStarFinder(fwhm = fwhm_mean + 0.5,minsep_fwhm = 3.0, threshold = 8*std,exclude_border = True,
                                         sharplo = 0.5,sharphi = 2.0, roundlo = 0.0,roundhi = 0.7)
            sources = daofind(data - bkg.background)

            #提取二次找星后使用的半径
            fwhm_cope = sources['fwhm']
            cope_mean,cope_median,cope_std = sigma_clipped_stats(fwhm_cope,sigma=3.0)


            radii = [cope_mean,cope_mean + 0.5]
            #判断星点数量是否够多
            if(len(sources) > star_num):
                positions = [(ix,iy) for ix,iy in zip(sources['xcentroid'],sources['ycentroid'])]
                #测光作业
                apertures = [pht.CircularAperture(positions,r = r * 2.5) for r in radii]
                aper_phot = pht.aperture_photometry(data - bkg.background,apertures, error = error)
                #将流量转化为星等
                for j in range(len(radii)):
                    fcol = 'aperture_sum_' + str(j)
                    ecol = 'aperture_sum_err_' + str(j)
                    flux = aper_phot[fcol]
                    fluxerr = aper_phot[ecol]
                    mag = -2.5 * np.log10(flux) + 25
                    magerr = 2.5 / (flux * np.log(10)) * fluxerr
                    aper_phot[fcol] = mag
                    aper_phot[ecol] = magerr
                    #重命名列名
                    aper_phot.rename_column(fcol,'mag_' + str(j))
                    aper_phot.rename_column(ecol,'magerr_' + str(j))

                #坐标转化
                x_pix = aper_phot['xcenter']
                y_pix = aper_phot['ycenter']
                ra_degree = w.pixel_to_world(x_pix, y_pix).ra.degree
                dec_degree = w.pixel_to_world(x_pix, y_pix).dec.degree

                #增加赤经赤纬列
                aper_phot['RA'] = ra_degree
                aper_phot['DEC'] = dec_degree
                #重命名列名
                #aper_phot.rename_column('xcenter','RA')
                #aper_phot.rename_column('ycenter','DEC')

                #为表单创建新的fits头
                new_header = fits.Header()
                new_header['EPOCH'] = raw_head['EPOCH']
                new_header['DATE-OBS'] = raw_head['DATE-OBS']                                                         
                if(bool_head == True):
                    new_header['TIME-OBS'] = raw_head['TIME-OBS']
                new_header['IMG_WIDTH'] = raw_head['NAXIS1']
                new_header['IMG_HEIGHT'] = raw_head['NAXIS2']
                bintable_hdu = fits.BinTableHDU(aper_phot,new_header)
                hdu1 = fits.HDUList([fits.PrimaryHDU(),bintable_hdu])
                hdu1.writeto(catfile,overwrite=True)
                #计数器+1
                progress = (i+1) / len(cfiles) * 100
                print(f"working process: [{int(progress):3d}%],",end='\r')
            else:
                #print('星点数量过少！疑似低质图像！' + rootname)
                error_image2.append(rootname)
                progress = (i+1) / len(cfiles) * 100
                print(f"working process: [{int(progress):3d}%],",end='\r')
                        
    
        print(" Completed   ! ")

    return error_image2