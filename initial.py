import os

# 子文件夹控制函数
def folder_control(reduced_folder):
    # 定义子文件夹名称
    subfolders = ['confirmed_star', 'png_images']

    # 遍历子文件夹列表查找是否存在子文件夹，不存在则创建
    for subfolder in subfolders:
        subfolder_path = os.path.join(reduced_folder, subfolder)  # 构建子文件夹路径
        if not os.path.exists(subfolder_path):  # 检查子文件夹是否存在
            os.makedirs(subfolder_path)  # 创建子文件夹
            #print(f"已创建子文件夹: {subfolder_path}")
        else:
            continue
            #print(f"子文件夹已存在: {subfolder_path}")



def address_cope (band_info,target_name):

    print("\033[1;31m Initialing file path!\033[0m")
    
    # 定义主文件夹路径
    reduced_folder = '../reduced_data/'
    #定义数据文件夹
    raw_folder = '../raw_data/'
    #调用文件夹控制函数
    folder_control(reduced_folder)


    #原始本底地址
    raw_bias = raw_folder +'/bias_*.fit'
    #超级本底地址
    super_bias = reduced_folder + 'superbias.fits'
    #改正后数据存储位置
    reduced_data_file_path = reduced_folder
    #图像存储地址
    png_file_path = reduced_folder + 'png_images/'
    #参考星与目标星存储地址
    confirmed_star_path = reduced_folder + 'confirmed_star/'

    #原始平场所在位置
    raw_flat = []
    #原始数据所在位置
    raw_data_file_path = []
    #超级平场所在位置
    super_flat = []
    #校正图像所在位置
    cor_data_file_path = []
    #裁剪校正图像所在位置
    cut_data_file_path = []
    #包含天测头的图像所在位置
    ast_data_file_path = []
    #星表地址
    star_cat_file_path = []
    #用于测试增益的本底地址
    biasfile1 = raw_folder + 'bias_001.fit'
    biasfile2 = raw_folder + 'bias_002.fit'
    #用于测试增益的平场地址
    flatfile1 = []
    flatfile2 = []
    #参考星所在目录
    confirmed_star_candidate_path = []
    #目标星所在目录
    target_star_path = []

    for i in range(len(band_info)):
        raw_flat.append(raw_folder + 'flat_' + band_info[i] + '_*.fit')
        super_flat.append(reduced_folder + 'superflat_' + band_info[i] + '.fits')
        raw_data_file_path.append(raw_folder+target_name+'_'+band_info[i]+'_*.fit')
        cor_data_file_path.append(reduced_folder+'Cor_'+target_name+'_'+band_info[i]+'_*.fits')
        cut_data_file_path.append(reduced_folder+'Cut_Cor_'+target_name+'_'+band_info[i]+'_???.fits')
        ast_data_file_path.append(reduced_folder+'Cut_Cor_'+target_name+'_'+band_info[i]+'_???-ast.fits')
        star_cat_file_path.append(reduced_folder+'Cut_Cor_'+target_name+'_'+band_info[i]+'_*-ast-cat.fits')
        flatfile1.append(raw_folder+'flat_'+band_info[i]+'_001.fit')
        flatfile2.append(raw_folder+'flat_'+band_info[i]+'_002.fit')
        confirmed_star_candidate_path.append(confirmed_star_path + band_info[i] + '_band_confirmed_star_*.fits')
        target_star_path.append(confirmed_star_path + band_info[i] + '_band_target_star.fits')

    return raw_bias,raw_data_file_path,super_bias,raw_flat,super_flat,cor_data_file_path,cut_data_file_path,ast_data_file_path,star_cat_file_path,confirmed_star_path,reduced_data_file_path,png_file_path,biasfile1,biasfile2,flatfile1,flatfile2,confirmed_star_candidate_path,target_star_path
