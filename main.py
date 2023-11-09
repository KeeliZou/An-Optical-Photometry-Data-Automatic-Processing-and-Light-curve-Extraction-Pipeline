#用户需要做的事：开辟一个文件夹，在里面创建raw_data,reduced_data,code三个文件夹
#把原始数据扔到raw_data里，把代码扔到reduced_data里

#这是需要用户自己输入的信息
#波段信息
band_info = ['B','V','R','I']
#目标名（与望远镜观测文件名中的名字一致）
target_name = 'V670And'
#以'07 35 29.88 +49 28 41.57'或'07:35:29.88 +49:28:41.57'的形式给出
target_coor = '23 30 37.2536 +46 24 04.2822'
#x、y轴裁剪值
edge_cutx = 20
edge_cuty = 20
#天测找星sigma倍数(默认50倍)
ast_sigma = 50
#星象初始半高全宽
fwhm_initial = 5
#图像星点阈值（低于此值意味着丢弃该图像）
star_num = 20
#绘图指定参考星索引（可以不设置）





#以下内容勿动
import initial
import bfc_2
import bfc_1
import bkg
import star_choose as sc
import plot as pt

#地址赋值
raw_bias,raw_data_file_path,super_bias,raw_flat,super_flat,cor_data_file_path,cut_data_file_path,ast_data_file_path,star_cat_file_path,confirmed_star_path,reduced_data_file_path,png_file_path,biasfile1,biasfile2,flatfile1,flatfile2,confirmed_star_candidate_path,target_star_path = initial.address_cope(band_info,target_name)

#本底合并
bfc_1.bias_combine(raw_bias,super_bias)
#平场合并
bfc_1.flat_combine(raw_flat,super_flat,super_bias,band_info)
#原始图像本底平场合并
bfc_2.bfc_raw(super_bias,super_flat,band_info,raw_data_file_path,reduced_data_file_path)
#图像裁剪
bfc_2.cor_cut(cor_data_file_path,reduced_data_file_path,band_info,edge_cutx,edge_cuty)
#评估天光背景，剔出不佳图像
cfiles = bkg.bkg_evaluate(band_info,cut_data_file_path)
#进行天测操作
error_image1 = bkg.ast_solve(band_info,cfiles,fwhm_initial,ast_sigma)
#计算增益与读出噪声
gain,rdnoise = bkg.gain_caculate(band_info,biasfile1,biasfile2,flatfile1,flatfile2)
#进行找星操作
error_image2 = bkg.star_find(band_info,ast_data_file_path,gain,fwhm_initial,star_num)
#这里看看要不要加一个error_image的保存
#星点对齐与参考星选择（外加望远镜赤经赤纬移动）
sta_image = sc.star_choose(target_coor,band_info,star_cat_file_path,confirmed_star_path,png_file_path)
#参考星选取
ref_star_cat,check_star_cat = sc.refer_chooes(band_info,confirmed_star_candidate_path,target_star_path,sta_image,png_file_path)
#进行绘图
pt.auto_plot(band_info,target_star_path,ref_star_cat,check_star_cat,confirmed_star_candidate_path,png_file_path,confirmed_star_path)