# -*- coding: utf-8 -*-
import xlrd
import xlwt 
import numpy as np

# 矩阵切片---Matrix slice(从属于ReadData的子函数)
def Matrix_slice(Matrix,Coordinate):
    row_start = Coordinate[0] #
    row_end   = Coordinate[1]
    col_start = Coordinate[2]
    col_end   = Coordinate[3]
    Matrix_partitioned = []  # A partitioned matrix
    for i in range(row_end-row_start):
        Matrix_partitioned.append([])
        for j in range(col_end-col_start):
            Matrix_partitioned[i].append(Matrix[row_start+i][col_start+j])
    return Matrix_partitioned

# ReadData---读取数据
def ReadData(filename):
    Data_origin = []
    readbook = xlrd.open_workbook(filename) # 打开excel
    sheetnum = len(readbook.sheets()) # 工作表数量
    # Data preprocessing
    for i in range(sheetnum):  # sheet number循环
        sheet = readbook.sheet_by_index(i) 
        n_row = sheet.nrows # 行数
        n_col = sheet.ncols # 列数
        Coordinate = [1,n_row,0,n_col]  # 坐标切片---coordinate of slice（要截取的范围，除了第一行，其他都要）
        Data_temp = sheet._cell_values  # excel表中的数据---data in the Excel file
        Data_origin.append(np.array(Matrix_slice(Data_temp,Coordinate)))
    return Data_origin


def WriteData(filename,x_node,ρ_A_b,ρ_B_b,ρ_C_b):
    workbook = xlwt.Workbook(encoding = 'utf-8')
    # 创建一个worksheet
    worksheet = workbook.add_sheet('x_node')
    # 写入excel
    for k in range(len(x_node)):
        worksheet.write(k+1,0,x_node[k])
    worksheet = workbook.add_sheet('y_b')
    for k in range(len(ρ_A_b)):
        worksheet.write(k+1,0,ρ_A_b[k])
        worksheet.write(k+1,1,ρ_B_b[k])
        worksheet.write(k+1,2,ρ_C_b[k])
    workbook.save(filename)