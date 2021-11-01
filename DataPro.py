import numpy as np
import copy

# %% 初始化数据
# 折现率计算函数
def Discount_Rate(r,LT):
    beta = r*(1 + r)**LT/((1 + r)**LT - 1)
    return beta

# 字符串是否含有数字判断函数
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    
class Parameter(object):
    def __init__(self,Data):
        # %% 其他参数
        self.r = 0.1 # 折现率
        self.Margin = 0.1 # 裕度
        
        self.Pr = 0.0922 # 主网购电成本（0.092k$/MWh）
        self.Outage_Pr = 1 # Outage cost（k$/MWh）
        self.C_ACDG = 0.0922 # ACDG   【WT】 k$/MWh
        self.S_ACDG = 0    # ACDG subsidy 
        self.C_DCDG = 0.159  # DCDG cost 【PV】k$/MWh
        self.S_DCDG = 0.106  # DCDG subsidy 
        
        # %% 可靠性管理参数
        self.SAIDI_ub = 10 # Reliability Regulation upper bound
        self.SAIDI_lb = 0.8 # Reliability Regulation lower bound
        self.SAIFI_ub = 1
        self.SAIFI_lb = 0.5
        
        # %% 潮流计算相关的基础参数
        self.Vr_AC = 4.16 # 交流额定电压（kV）
        self.Vr_DC = 4.16/(3 / 8)**0.5 # 直流额定电压（kV）
        self.AClineSmax = 10 # 交流线路容量S（MVA）
        self.DClineSmax = 10 # 直流线路容量S（MVA）
        self.Vmax_AC = 1.05*self.Vr_AC # 交流最大电压（kV）
        self.Vmin_AC = 0.85*self.Vr_AC # 交流最小电压（kV） # 修改测试说明！！
        self.Vmax_DC = 1.05*self.Vr_DC # 直流最大电压（kV）
        self.Vmin_DC = 0.85*self.Vr_DC # 直流最小电压（kV） # 修改测试说明！！
        self.Imax_AC =  self.AClineSmax / self.Vr_AC / 3**0.5 # 交流最大电流（kA）
        self.Imax_DC =  self.DClineSmax / self.Vr_DC # 直流最大电流（kA）
        
        # %% 设备基础参数
        self.η_INV = 0.98 # 逆变器效率
        self.η_REC = 0.98 # 整流器效率
        self.η_VSC = 0.98 # VSC效率
        self.VSC_Kc = (3 / 8)**0.5 # VSC的Kc值
        # self.VSC_Mc = 0.96 # VSC的Mc值
        self.VSC_Mc = 0.99 # VSC的Mc值
        self.INV_Kc = (3 / 8)**0.5 # INV的Kc值
        # self.INV_Mc = 0.96 # INV的Mc值
        self.INV_Mc = 0.99 # INV的Mc值
        self.REC_Kc = (3 / 8)**0.5 # REC的Kc值
        # self.REC_Mc = 0.96 # REC的Mc值
        self.REC_Mc = 0.99 # REC的Mc值
        
        # %% 设备生命周期
        self.LT_VSC = 20 # VSC LifeTime
        self.LT_REC = 20 # Rectifier LifeTime
        self.LT_INV = 20 # Inverter LifeTime
        self.LT_ACB = 20 # AC breaker LifeTime
        self.LT_DCB = 20 # DC breaker LifeTime
        self.LT_ACL = 20 # AC Line LifeTime
        self.LT_DCL = 20 # DC Line LifeTime
        
        
        # %% 设备成本参数
        self.c_VSC = 0.170 # 170k$/MVA
        self.c_INV = 0.170 # 170k$/MVA
        self.c_REC = 0.170 # 170k$/MVA
        self.C_ACB = 24    # 24K$
        self.C_DCB = 26    # 26K$ 
        # self.C_DCB = 2.6    # 24K$
        self.c_ACL = 17.5  # AC Line k$/km
        self.c_DCL = 11.6  # DC Line k$/km
        # self.c_DCL = 1.16  # DC Line
        
        # %% 设备维修成本参数
        self.OC_Line_A = 0.45 # Line A 运维成本
        self.OC_Line_B = 0.3  # Line B 运维成本
        self.OC_Line_C = 0.3  # Line C 运维成本
        

        # %% 大表格数据
        self.Branch = Data[0] # 支路数据
        self.Node = Data[1] # 节点数据
        
        # %% 环路数据
        tmp = Data[2] # 环路数据
        self.LoopNum = np.size(tmp,0)
        col = np.size(tmp,1)
        Loop = [[] for i in range(self.LoopNum)]
        for i in range(self.LoopNum):
            for j in range(col):
                if is_number(tmp[i,j]):
                    Loop[i].append(int(float(tmp[i,j])))
        self.Loop = Loop # 处理为二维list
         
        # %% 潮流计算节点参数
        self.AC_Demand = self.Node[:,1] # 交流负荷有功峰值/MW
        self.DC_Demand = self.Node[:,2] # 交流负荷有功峰值/MW
        self.AC_DG_Cap = self.Node[:,3] # 交流DG容量/MW
        self.DC_DG_Cap = self.Node[:,4] # 直流DG容量/MW
        
        

        # %% 潮流计算支路参数
        self.L_Branch = self.Branch[:,3] # 支路长度

        
        # %% 常用数据
        self.N_Node = np.size(self.Node,0) # 节点数
        self.N_Branch = np.size(self.Branch,0) # 支路数
        self.N_Customer = self.Node[:,5] # 节点用户数
        self.N_AC_Customer = self.Node[:,6] # 交流节点负荷数
        self.N_DC_Customer = self.Node[:,7] # 直流流节点负荷数

        
        
        # %% AD/DC线路阻抗参数
        self.ACL_Z = 0.5287 # AC线路单位阻抗
        self.DCL_Z = 0.2744 # DC线路单位阻抗
        # self.DCL_Z = 0 # DC线路单位阻抗
        self.Z_Branch_A = np.zeros(self.N_Branch)
        self.Z_Branch_B = np.zeros(self.N_Branch)
        self.Z_Branch_C = np.zeros(self.N_Branch)
        for k in range(self.N_Branch):
            self.Z_Branch_A[k] = self.L_Branch[k] * self.ACL_Z
            self.Z_Branch_B[k] = self.L_Branch[k] * self.DCL_Z
        self.Z_Branch_C = self.Z_Branch_B    
            
        # %% 节点-支路关联矩阵
        A_Matrix = np.zeros((self.N_Node, self.N_Branch)) # 节点数*支路数
        for k in range(self.N_Branch):
            A_Matrix[int(self.Branch[k,1]),k] = 1
            A_Matrix[int(self.Branch[k,2]),k] = -1
        self.A_Matrix = A_Matrix 
        
        # %% 元件可靠性参数
        K_c = 1 # 可分析
        K_s = 1 # 可分析
        self.DCL_FaultRate = np.zeros(self.N_Branch) # 单位长度直流线路故障率
        self.DCL_RepairTime = np.zeros(self.N_Branch) # 单位长度直流线路故障修复时间
        self.ACL_FaultRate = np.zeros(self.N_Branch) # 单位长度交流线路故障率
        self.ACL_RepairTime = np.zeros(self.N_Branch) # 单位长度交流线路故障修复时间
        self.dcl_faultrate = 0.335107/1.609 # 直流线路单位长度的故障率
        # self.dcl_repairtime = 4 # 直流线路单位长度的修复时间 
        self.dcl_repairtime = 4/1.609 # 直流线路单位长度的修复时间 
        self.acl_faultrate = 0.251330/1.609 # 交流线路单位长度的故障率
        # self.acl_repairtime = 4 # 交流线路单位长度的修复时间 
        self.acl_repairtime = 4/1.609 # 直流线路单位长度的修复时间 
        
        self.VSC_FaultRate = 1.4 * K_s # VSC故障率
        self.IC_FaultRate = 0.660515 * K_s # 逆变器故障率
        self.RC_FaultRate = 2.591962 * K_s # 整流器故障率
        
        self.VSC_RepairTime = 4.1 # VSC故障修复时间 
        self.IC_RepairTime = 5.321 # 逆变器故障修复时间
        self.RC_RepairTime = 3.471 # 整流器故障修复时间 
        
        self.Du_A_b = np.zeros(self.N_Branch) # A类线路故障时的等效停电时间
        self.Du_B_b = np.zeros(self.N_Branch) # B类线路故障时的等效停电时间
        self.Du_C_b = np.zeros(self.N_Branch) # C类线路故障时的等效停电时间
        
        # 断路器
        self.ACB_FaultRate = 0.006 * K_c
        self.DCB_FaultRate = 0.2999 * K_c
        self.ACB_RepairTime = 4 * K_c
        self.DCB_RepairTime = 10 * K_c
        
        
        # %% 可靠性参数运算
        for k in range(self.N_Branch):
            self.DCL_FaultRate[k] = self.dcl_faultrate * self.Branch[k,3] # 直流线路故障率
            self.DCL_RepairTime[k] = self.dcl_repairtime * 1 # 直流线路故障修复时间
            self.ACL_FaultRate[k] = self.acl_faultrate * self.Branch[k,3] # 交流线路故障率
            self.ACL_RepairTime[k] = self.acl_repairtime * 1 # 交流线路故障修复时间
            self.Du_A_b[k] = self.ACL_FaultRate[k] * self.ACL_RepairTime[k] + 2*self.ACB_FaultRate # A型线路等效故障时间
            self.Du_B_b[k] = self.DCL_FaultRate[k] * self.DCL_RepairTime[k] + self.VSC_FaultRate * self.VSC_RepairTime  + self.ACB_FaultRate * self.ACB_RepairTime + 2*self.DCB_FaultRate * self.DCB_RepairTime# B型线路等效故障时间
            self.Du_C_b[k] = self.DCL_FaultRate[k] * self.DCL_RepairTime[k] + 2*self.DCB_FaultRate # C型线路等效故障时间
        
        # %% 节点分类
        self.S_Node = self.Node[self.Node[:,8] == 1,0].astype(int) # 电源节点
        self.L_Node = self.Node[self.Node[:,8] == 0,0].astype(int) # 负荷节点
        self.N_S_Node = len(self.S_Node) # 电源点数量
        self.Node2Substation = self.Node[:,12] # 节点归属于哪个变电站
        
        
        # %% 负荷和DG的存在系数
        self.exist_AC_Demand = self.Node[:,1] != 0 # AC负荷的存在系数
        self.exist_DC_Demand = self.Node[:,2] != 0 # DC负荷的存在系数
        self.exist_AC_DG_Cap = self.Node[:,3] != 0 # AC负荷的存在系数
        self.exist_DC_DG_Cap = self.Node[:,4] != 0 # DC负荷的存在系数
        
        # %% 每个变电站管辖的节点计算
        self.GN_s= np.zeros(self.N_Node)
        self.GDD_s= np.zeros(self.N_Node)
        self.GDDGC_s= np.zeros(self.N_Node)
        for s in self.S_Node:
            self.GN_s[s] = sum(self.N_Customer[self.Node2Substation[:] == s])  # 每个电源下的总用户数
            self.GDD_s[s] = sum(self.AC_Demand[self.Node2Substation[:] == s]) + sum(self.DC_Demand[self.Node2Substation[:] == s])   # 每个电源下的总负荷
            self.GDDGC_s[s] = sum(self.AC_DG_Cap[self.Node2Substation[:] == s]) + sum(self.DC_DG_Cap[self.Node2Substation[:] == s])  # 每个电源下的总负荷
        
        self.GN = sum(self.GN_s) # 总用户数
        self.GDD = sum(self.GDD_s) # 总负荷
        self.GDDGC = sum(self.AC_DG_Cap) + sum(self.DC_DG_Cap) # 总DG
        self.AP = self.GDD/sum(self.N_Customer) # 平均负荷
        
        # %% 画图坐标数据
        self.Node_x = self.Node[:,10] # 节点x坐标
        self.Node_y = self.Node[:,11] # 节点y坐标
        
        # %% 建设变量 （已知参数）
        self.x_node_type = self.Node[:,9]
        self.ρ_A_b = self.Branch[:,4]
        self.ρ_B_b = self.Branch[:,5]
        self.ρ_C_b = self.Branch[:,6]
        
        
        # %% 成本参数运算
        self.IC_Bus_AC_i = np.zeros(self.N_Node)
        self.IC_Bus_DC_i = np.zeros(self.N_Node)
        
        self.IC_Line_A = np.zeros(self.N_Branch)
        self.IC_Line_B_p1 = np.zeros(self.N_Branch)
        self.IC_Line_B_p2 = np.zeros(self.N_Branch)
        self.IC_Line_C = np.zeros(self.N_Branch)
        
        #  IC_Bus_AC
        for i in range(self.N_Node):
            self.IC_Bus_AC_i[i] = Discount_Rate(self.r,self.LT_REC) * self.c_REC * (1 + self.Margin) * self.DC_Demand[i] \
            + Discount_Rate(self.r,self.LT_INV) * self.c_INV * (1 + self.Margin) * self.DC_DG_Cap[i] \
            + self.exist_AC_Demand[i] * Discount_Rate(self.r,self.LT_ACB) * self.C_ACB \
            + self.exist_DC_Demand[i] * (Discount_Rate(self.r,self.LT_ACB) * self.C_ACB + Discount_Rate(self.r,self.LT_DCB) * self.C_DCB) \
            + self.exist_AC_DG_Cap[i] * Discount_Rate(self.r,self.LT_ACB) * self.C_ACB \
            + self.exist_DC_DG_Cap[i] * (Discount_Rate(self.r,self.LT_ACB) * self.C_ACB + Discount_Rate(self.r,self.LT_DCB) * self.C_DCB)
        #  IC_Bus_DC
        for i in range(self.N_Node):
            self.IC_Bus_DC_i[i] = Discount_Rate(self.r,self.LT_INV) * self.c_INV * (1 + self.Margin) * self.AC_Demand[i] \
            + Discount_Rate(self.r,self.LT_REC) * self.c_REC * (1 + self.Margin) * self.AC_DG_Cap[i] \
            + self.exist_AC_Demand[i] * (Discount_Rate(self.r,self.LT_ACB) * self.C_ACB + Discount_Rate(self.r,self.LT_DCB) * self.C_DCB) \
            + self.exist_DC_Demand[i] * Discount_Rate(self.r,self.LT_DCB) * self.C_DCB \
            + self.exist_AC_DG_Cap[i] * (Discount_Rate(self.r,self.LT_ACB) * self.C_ACB + Discount_Rate(self.r,self.LT_DCB) * self.C_DCB) \
            + self.exist_DC_DG_Cap[i] * Discount_Rate(self.r,self.LT_DCB) * self.C_DCB
        #  IC_Line_A
        for b in range(self.N_Branch):
            self.IC_Line_A[b] = Discount_Rate(self.r,self.LT_ACL) * self.c_ACL * self.L_Branch[b] + 2 * Discount_Rate(self.r,self.LT_ACL) * self.C_ACB
            
        
        #  IC_Line_B
        for b in range(self.N_Branch):
            self.IC_Line_B_p1[b] = ( Discount_Rate(self.r,self.LT_DCL) * self.c_DCL * self.L_Branch[b] + Discount_Rate(self.r,self.LT_ACB) * self.C_ACB
            + 2 * Discount_Rate(self.r,self.LT_DCB) * self.C_DCB )
    
        for b in range(self.N_Branch):
            self.IC_Line_B_p2[b] = ( Discount_Rate(self.r,self.LT_VSC) * self.c_VSC * (1 + self.Margin) * self.Vr_DC )
        
        #  IC_Line_C
        for b in range(self.N_Branch):
            self.IC_Line_C[b] = Discount_Rate(self.r,self.LT_DCL) * self.c_DCL * self.L_Branch[b] + 2 * Discount_Rate(self.r,self.LT_DCB) * self.C_DCB

        # %% 变量的初始索引(定义的时候也要按顺序定义！)
        self.V_index = 0
        self.I_index = self.V_index + self.N_Node
        self.x_type_i_index = self.I_index + self.N_Branch
        self.y_b_index = self.x_type_i_index + self.N_Node
        self.ρ_A_b_index = self.y_b_index + self.N_Branch
        self.ρ_B_b_index = self.ρ_A_b_index + self.N_Branch
        self.ρ_C_b_index = self.ρ_B_b_index + self.N_Branch
        self.SAIFI_index = self.ρ_C_b_index + self.N_Branch
        self.SAIDI_index = self.SAIFI_index + 1
        
        # %% 网损计算相关
        self.c_L = 0.04*8760 # 网损系数 0.04k/MWh
        self.N_v = 4 # 分段数,这尼玛扯淡啊
        self.Imax = max(self.Imax_AC,self.Imax_DC)
        # self.Imax = 1
        self.d_v = 2 * self.Imax / self.N_v
        self.I_b_v = [- self.Imax + v * self.d_v for v in range(self.N_v+1)]
        self.s_v = np.zeros(self.N_v)
        for v in range(self.N_v):
            self.s_v[v] = self.I_b_v[v+1] + self.I_b_v[v]