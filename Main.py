'''
Developer: Ziyao Wang
Place: College of Electric Power, South China University of Technology
Content: Reliability-based AC-DC Distritbution Network Planning Model
Method: MILP
Time: 2021/9/13
'''
from gurobipy import *
import numpy as np
import DataGen
import DataPro
import PlotFunc
from pandas import Series,DataFrame
import pandas as pd



if __name__ == '__main__':
    # filename = "Data_6 - LowDG.xlsx" # 1-3
    filename = "Data_6 - HighDG.xlsx" # 4
    Data = DataGen.ReadData(filename)
    Para = DataPro.Parameter(Data)
    
    model = Model()
    
    # %% Define Variable
    r = 0.1 # DiscountRate
    Margin = 0.1 # Margin 
    km = 0.1 # coefficient 
    # Line Type Variable
    ρ_A_b = model.addVars(Para.N_Branch,vtype = GRB.BINARY, name = 'rou_A_b') # equal to 1 when the branch b is A Line, being 0 otherwise
    ρ_B_b = model.addVars(Para.N_Branch,vtype = GRB.BINARY, name = 'rou_B_b')
    ρ_C_b = model.addVars(Para.N_Branch,vtype = GRB.BINARY, name = 'rou_C_b')
    
    # Node Type Variable
    x_type_i = model.addVars(Para.N_Node,vtype = GRB.BINARY, name = 'x_type_i') # equal to 1 when the bus is DC bus, being 0 when the bus is AC bus
    
    # Line Investment Variable
    y_b = model.addVars(Para.N_Branch,vtype = GRB.BINARY, name = 'y_b') # equal to 1 when the branch is built, being 0 otherwise
                                                                                                                                                                                                                                                                                                                                                                                                                                     
    # Intermediate Variable for Line
    ρ_A_b_tmp = model.addVars(Para.N_Branch,vtype = GRB.BINARY, name = 'ρ_A_b_tmp')
    ρ_B_b_tmp = model.addVars(Para.N_Branch,vtype = GRB.BINARY, name = 'ρ_B_b_tmp') 
    ρ_C_b_tmp = model.addVars(Para.N_Branch,vtype = GRB.BINARY, name = 'ρ_C_b_tmp') 
    
    # Reliability Variable
    EENS = model.addVar(name = 'EENS') 
    BEENS_A = model.addVars(Para.N_Branch,name = 'BEENS_A') # BEENS of Type A Branch 
    BEENS_B = model.addVars(Para.N_Branch,name = 'BEENS_B') # BEENS of Type B Branch 
    BEENS_C = model.addVars(Para.N_Branch,name = 'BEENS_C') # BEENS of Type C Branch 
    NEENS = model.addVars(Para.N_Node, name = 'NEENS')

    SAIFI = model.addVar(name = 'SAIFI')
    BSAIFI_A = model.addVars(Para.N_Branch, name = 'BSAIFI_A')
    BSAIFI_B = model.addVars(Para.N_Branch, name = 'BSAIFI_B')
    BSAIFI_C = model.addVars(Para.N_Branch, name = 'BSAIFI_C')
    NSAIFI = model.addVars(Para.N_Node, name = 'NSAIFI')
    
    
    SAIDI = model.addVar(name = 'SAIDI')
    BSAIDI_A = model.addVars(Para.N_Branch, name = 'BSAIDI_A')
    BSAIDI_B = model.addVars(Para.N_Branch, name = 'BSAIDI_B')
    BSAIDI_C = model.addVars(Para.N_Branch, name = 'BSAIDI_C')
    NSAIDI = model.addVars(Para.N_Node, name = 'N_SAIDI')
    
    # CP
    CP_b = model.addVars(Para.N_Branch,lb = - GRB.INFINITY, name = 'CP_b') # CP_b: Curtailed power due to outage of branch b 
    α_CP = model.addVars(Para.N_Branch,vtype = GRB.BINARY, name = 'α_CP') # α_CP: Auxiliary Variable 
    
    # DD
    DD_b = model.addVars(Para.N_Branch,lb = - GRB.INFINITY, name = 'DD_b') # DD_b
    DD_b_p = model.addVars(Para.N_Branch,name = 'DD_b_p') # DD_b_p
    DD_b_n = model.addVars(Para.N_Branch,name = 'DD_b_n') # DD_b_n
    α_DD = model.addVars(Para.N_Branch,vtype = GRB.BINARY, name = 'α_DD ') # α_DD 
    
    # DDGC
    DDGC_b = model.addVars(Para.N_Branch,lb = - GRB.INFINITY, name = 'DDGC_b') # DDGC_b
    DDGC_b_p = model.addVars(Para.N_Branch,name = 'DDGC_b_p') # DDGC_b_p
    DDGC_b_n = model.addVars(Para.N_Branch,name = 'DDGC_b_n') # DDGC_b_n
    α_DDGC = model.addVars(Para.N_Branch,vtype = GRB.BINARY, name = 'α_DDGC ') # α_DDGC 
    
    # CCN
    CCN_b = model.addVars(Para.N_Branch,lb = - GRB.INFINITY,  name = 'CCN')
    α_CCN = model.addVars(Para.N_Branch, vtype = GRB.BINARY, name = 'α_CCN')
    
    # DCN
    DCN_b = model.addVars(Para.N_Branch, name = 'DCN_b')
    DCN_b_p = model.addVars(Para.N_Branch, name = 'DCN_b_p')
    DCN_b_n = model.addVars(Para.N_Branch, name = 'DCN_b_n')
    α_DCN = model.addVars(Para.N_Branch,vtype = GRB.BINARY, name = 'α_DCN')
    
    # Investment Cos
    COST = model.addVar(name = 'COST')
    C_I = model.addVar(name = 'C_I')
    C_M = model.addVar(name = 'C_M')
    C_R = model.addVar(name = 'C_R')
    C_E = model.addVar(name = 'C_E')
    
    C_EENS = model.addVar(name = 'C_EENS')
    C_RR = model.addVar(name = 'C_RR')
    C_RR_SAIDI = model.addVar(name = 'C_RR_SAIDI')
    C_RR_SAIFI = model.addVar(name = 'C_RR_SAIFI')
    
    # Investment Cost
    IC_Bus_AC_i = model.addVars(Para.N_Node, name = 'IC_Bus_AC')
    IC_Bus_DC_i = model.addVars(Para.N_Node, name = 'IC_Bus_DC')
    IC_Line_A = model.addVars(Para.N_Branch, name = 'IC_Line_A')
    IC_Line_B_p1 = model.addVars(Para.N_Branch, name = 'IC_Line_B_p1')
    IC_Line_B_p2 = model.addVars(Para.N_Branch, name = 'IC_Line_B_p2')
    IC_Line_C = model.addVars(Para.N_Branch, name = 'IC_Line_C')
    
    
    # Powerflow Variable
    I_b = model.addVars(Para.N_Branch,lb = -GRB.INFINITY, name = 'I_b')
    V_i = model.addVars(Para.N_Node, name = 'V_i')
    δ_i = model.addVars(Para.N_Node,lb = - GRB.INFINITY,  name = 'δ_i')
    v_b = model.addVars(Para.N_Branch,lb = -GRB.INFINITY, name = 'v_b')
    P_s = model.addVars(Para.N_Node,lb = - GRB.INFINITY,  name = 'P_s')
    
    # Intermedian Variable
    u1 = model.addVars(4,  name = 'u1')
    u2 = model.addVars(4,  name = 'u2')
    sita1 = model.addVars(3,vtype = GRB.BINARY, name = 'sita1')
    sita2 = model.addVars(3,vtype = GRB.BINARY, name = 'sita2')
    
    # Power Loss
    C_L = model.addVar(name = 'C_L')
    deta_b_v = model.addVars(Para.N_Branch, Para.N_v, ub = Para.d_v, name = 'deta_b_v')
    deta_b_v_A = model.addVars(Para.N_Branch, Para.N_v, name = 'deta_b_v_A')
    deta_b_v_B = model.addVars(Para.N_Branch, Para.N_v, name = 'deta_b_v_B')
    deta_b_v_C = model.addVars(Para.N_Branch, Para.N_v, name = 'deta_b_v_C')
    
    # %% Objective function: Investment Cost + Maintenace Cost + Eletricity Production Cost + Reliability Cost
    model.addConstr( COST == C_I + C_M + C_E + C_R) # (1)
    
    C_Eco = model.addVar(name = 'C_Eco')
    model.addConstr( C_Eco == C_I + C_M + C_E) 
    C_Rel = model.addVar(name = 'C_Rel')
    model.addConstr( C_Rel == C_R ) 
    
    
    # %% (1) Investment Cost
    model.addConstr(C_I == quicksum(Para.IC_Bus_AC_i[i] * (1-x_type_i[i]) for i in range (Para.N_Node))
                    + quicksum(Para.IC_Bus_DC_i[i] * x_type_i[i] for i in range (Para.N_Node))
                    + quicksum(Para.IC_Line_A[b] * ρ_A_b[b] for b in range(Para.N_Branch))
                    + quicksum(Para.IC_Line_C[b] * ρ_C_b[b] for b in range(Para.N_Branch))
                    + quicksum((Para.IC_Line_B_p1[b] * ρ_B_b[b] + Para.IC_Line_B_p2[b] * v_b[b]) for b in range(Para.N_Branch))) # !!!
    
    # %% (2) Maintenance Cost
    model.addConstr(
        C_M == quicksum((Para.OC_Line_A * ρ_A_b[b] + Para.OC_Line_B * ρ_B_b[b] + Para.OC_Line_C * ρ_C_b[b]) for b in range(Para.N_Branch))
        + quicksum(Para.IC_Bus_AC_i[i] * km * (1-x_type_i[i]) for i in range(Para.N_Node))
        + quicksum(Para.IC_Bus_DC_i[i] * km * x_type_i[i] for i in range(Para.N_Node))
        )
    
    # %% (3) Electricity Production Cost
    model.addConstr(
        C_E == quicksum( 8760 * Para.Pr * P_s[s] for s in Para.S_Node ) 
        + quicksum( 8760 * (Para.C_ACDG - Para.S_ACDG) * Para.AC_DG_Cap[i]  for i in range(Para.N_Node))
        + quicksum( 8760 * (Para.C_DCDG - Para.S_DCDG) * Para.DC_DG_Cap[i]  for i in range(Para.N_Node))
        )
    
    # %% (4) Reliability Cost
    # ks = 1
    # kc = 0.9
    model.addConstr(C_R ==  C_EENS + C_RR)
    model.addConstr(C_EENS == EENS * Para.Outage_Pr )
    model.addConstr(C_RR == C_RR_SAIDI + C_RR_SAIFI )

    model.addConstr(C_RR_SAIDI == 0.3 * SAIDI )
    model.addConstr(C_RR_SAIFI == 1  * SAIFI )
    
    
    # %% (5) Power Loss Cost
    # miu_A = model.addVars(Para.N_Branch, Para.N_v, name = 'miu_A')
    # miu_B = model.addVars(Para.N_Branch, Para.N_v, name = 'miu_B')
    # miu_C = model.addVars(Para.N_Branch, Para.N_v, name = 'miu_C')
    
    # M_Loss = 5
    # for b in range(Para.N_Branch):
    #     for v in range(Para.N_v):
    #         model.addConstr(miu_A[b,v] >= -M_Loss * ρ_A_b[b])
    #         model.addConstr(miu_A[b,v] <= M_Loss * ρ_A_b[b])
    #         model.addConstr(miu_A[b,v] >= -M_Loss * (1 - ρ_A_b[b]) + deta_b_v[b,v])
    #         model.addConstr(miu_A[b,v] <= M_Loss * (1 - ρ_A_b[b]) + deta_b_v[b,v])
    
    # for b in range(Para.N_Branch):
    #     for v in range(Para.N_v):
    #         model.addConstr(miu_B[b,v] >= -M_Loss * ρ_B_b[b])
    #         model.addConstr(miu_B[b,v] <= M_Loss * ρ_B_b[b])
    #         model.addConstr(miu_B[b,v] >= -M_Loss * (1 - ρ_B_b[b]) + deta_b_v[b,v])
    #         model.addConstr(miu_B[b,v] <= M_Loss * (1 - ρ_B_b[b]) + deta_b_v[b,v])
    
    # for b in range(Para.N_Branch):
    #     for v in range(Para.N_v):
    #         model.addConstr(miu_C[b,v] >= -M_Loss * ρ_C_b[b])
    #         model.addConstr(miu_C[b,v] <= M_Loss * ρ_C_b[b])
    #         model.addConstr(miu_C[b,v] >= -M_Loss * (1 - ρ_C_b[b]) + deta_b_v[b,v])
    #         model.addConstr(miu_C[b,v] <= M_Loss * (1 - ρ_C_b[b]) + deta_b_v[b,v])
            

    # I_2 = model.addVars(Para.N_Branch, name = 'I_2')
    # PL_A = model.addVars(Para.N_Branch, name = 'PL_A')
    # PL_B = model.addVars(Para.N_Branch, name = 'PL_B')
    # PL_C = model.addVars(Para.N_Branch, name = 'PL_C')

    # for b in range(Para.N_Branch):
    #     model.addConstr(I_b[b] == -Para.Imax + quicksum(deta_b_v[b,v] for v in range(Para.N_v)))
    #     model.addConstr(I_2[b] == Para.Imax**2 + quicksum(Para.s_v[i]*deta_b_v[b,i] for i in range(Para.N_v)))
    #     model.addConstr(PL_A[b] == Para.Z_Branch_A[b] * ρ_A_b[b] * Para.Imax**2 + Para.Z_Branch_A[b] * quicksum(Para.s_v[i]*miu_A[b,i] for i in range(Para.N_v)))
    #     model.addConstr(PL_B[b] == Para.Z_Branch_B[b] * ρ_B_b[b] * Para.Imax**2 + Para.Z_Branch_B[b] * quicksum(Para.s_v[i]*miu_B[b,i] for i in range(Para.N_v)))
    #     model.addConstr(PL_C[b] == Para.Z_Branch_C[b] * ρ_C_b[b] * Para.Imax**2 + Para.Z_Branch_C[b] * quicksum(Para.s_v[i]*miu_C[b,i] for i in range(Para.N_v)))
        
    # # model.addConstr(C_L == I_2[0])
    # model.addConstr(C_L == Para.c_L * quicksum( PL_A[b] + PL_B[b] + PL_C[b] for b in range(Para.N_Branch)) )
    
    # PL = model.addVar( name = 'PL')
    # model.addConstr(PL == quicksum( PL_A[b] + PL_B[b] + PL_C[b] for b in range(Para.N_Branch)) )
    
    '''
    I2_A = model.addVar(name = 'I2_A')
    I2_B = model.addVar(name = 'I2_B')
    I2_C = model.addVar(name = 'I2_C')
    
    model.addConstr(I2_A == quicksum(Para.Z_Branch_A[b] * ρ_A_b[b] * Para.Imax**2 + Para.Z_Branch_A[b] * quicksum(Para.s_v[v]*miu_A[b,v] for v in range(Para.N_v)) for b in range(Para.N_Branch)))
    model.addConstr(I2_B == quicksum(Para.Z_Branch_B[b] * ρ_B_b[b] * Para.Imax**2 + Para.Z_Branch_B[b] * quicksum(Para.s_v[v]*miu_B[b,v] for v in range(Para.N_v)) for b in range(Para.N_Branch)))
    model.addConstr(I2_C == quicksum(Para.Z_Branch_C[b] * ρ_C_b[b] * Para.Imax**2 + Para.Z_Branch_C[b] * quicksum(Para.s_v[v]*miu_C[b,v] for v in range(Para.N_v)) for b in range(Para.N_Branch)))

    model.addConstr(C_L == Para.c_L * ( 
                    quicksum(Para.Z_Branch_A[b] * ρ_A_b[b] * Para.Imax**2 + Para.Z_Branch_A[b] * quicksum(Para.s_v[v]*miu_A[b,v] for v in range(Para.N_v)) for b in range(Para.N_Branch))
                    + quicksum(Para.Z_Branch_B[b] * ρ_B_b[b] * Para.Imax**2 + Para.Z_Branch_B[b] * quicksum(Para.s_v[v]*miu_B[b,v] for v in range(Para.N_v)) for b in range(Para.N_Branch))
                    + quicksum(Para.Z_Branch_C[b] * ρ_C_b[b] * Para.Imax**2 + Para.Z_Branch_C[b] * quicksum(Para.s_v[v]*miu_C[b,v] for v in range(Para.N_v)) for b in range(Para.N_Branch))    
                    )
                    )
    '''
    
    # %% Constraint: Investment variable and investment type variable
    for k in range(Para.N_Branch):
        model.addConstr(y_b[k] == ρ_A_b[k] + ρ_B_b[k] + ρ_C_b[k])
        model.addConstr(y_b[k]  <= 1)
        model.addConstr(ρ_A_b[k]>=0)
        model.addConstr(ρ_A_b[k]<=ρ_A_b_tmp[k])
        model.addConstr(ρ_B_b[k]>=0)
        model.addConstr(ρ_B_b[k]<=ρ_B_b_tmp[k])
        model.addConstr(ρ_C_b[k]>=0)
        model.addConstr(ρ_C_b[k]<=ρ_C_b_tmp[k])
        
    for k in range(Para.N_Branch):
        i = int(Para.Branch[k,1])
        j = int(Para.Branch[k,2])
        model.addConstr(ρ_A_b_tmp[k] <= 1 - x_type_i[i]) 
        model.addConstr(ρ_A_b_tmp[k] <= 1 - x_type_i[j]) 
        model.addConstr(ρ_A_b_tmp[k] >= 1 - (x_type_i[i] + x_type_i[j])) 
        
        model.addConstr(ρ_B_b_tmp[k] <= (x_type_i[i] + x_type_i[j])) 
        model.addConstr(ρ_B_b_tmp[k] <= 2 - (x_type_i[i] + x_type_i[j])) 
        model.addConstr(ρ_B_b_tmp[k] >= x_type_i[i] - x_type_i[j])
        model.addConstr(ρ_B_b_tmp[k] >= x_type_i[j] - x_type_i[i]) 
        
        model.addConstr(ρ_C_b_tmp[k] <= x_type_i[i]) 
        model.addConstr(ρ_C_b_tmp[k] <= x_type_i[j]) 
        model.addConstr(ρ_C_b_tmp[k] >= (x_type_i[i] + x_type_i[j]) - 1) 
        
    # %% Constraint: Reliability Constraint
    # %% (1) EENS
    model.addConstr(EENS == quicksum(BEENS_A[k] for k in range(Para.N_Branch))
                      + quicksum(BEENS_B[k] for k in range(Para.N_Branch)) 
                      + quicksum(BEENS_C[k] for k in range(Para.N_Branch)) 
                      + quicksum(NEENS[i] for i in range(Para.N_Node)))
    
    # BEENS
    M_BEENS = 10*(Para.GDD)
    for k in range(Para.N_Branch):
        model.addConstr(BEENS_A[k] - Para.Du_A_b[k] * CP_b[k] <= M_BEENS * (1 - ρ_A_b[k]) ) 
        model.addConstr(BEENS_A[k] - Para.Du_A_b[k] * CP_b[k] >= -M_BEENS * (1 - ρ_A_b[k]) ) 
        model.addConstr(BEENS_A[k]  <= M_BEENS * ρ_A_b[k] ) 
        
        model.addConstr(BEENS_B[k] - Para.Du_B_b[k] * CP_b[k] <= M_BEENS * (1 - ρ_B_b[k]) ) 
        model.addConstr(BEENS_B[k] - Para.Du_B_b[k] * CP_b[k] >= -M_BEENS * (1 - ρ_B_b[k]) ) 
        model.addConstr(BEENS_B[k]  <= M_BEENS * ρ_B_b[k] ) 
        
        model.addConstr(BEENS_C[k] - Para.Du_C_b[k] * CP_b[k] <= M_BEENS * (1 - ρ_C_b[k]) ) 
        model.addConstr(BEENS_C[k] - Para.Du_C_b[k] * CP_b[k] >= -M_BEENS * (1 - ρ_C_b[k]) ) 
        model.addConstr(BEENS_C[k]  <= M_BEENS * ρ_C_b[k] ) 
    
    M_CP = Para.GDD * 2
    for k in range(Para.N_Branch):
        model.addConstr(CP_b[k] >= DD_b[k] - DDGC_b[k]) 
        model.addConstr(CP_b[k] <= DD_b[k] - DDGC_b[k] + M_CP *(1 - α_CP[k])) 
        model.addConstr(CP_b[k] >= 0) 
        model.addConstr(CP_b[k] <= M_CP *α_CP[k]) 
    
    # DD_b
    M_DD = Para.GDD * 2
    for k in range(Para.N_Branch):
        model.addConstr(DD_b[k] == DD_b_p[k] + DD_b_n[k]) 
        model.addConstr(DD_b_p[k] >= 0) 
        model.addConstr(DD_b_p[k] <= M_DD*α_DD[k]) 
        model.addConstr(DD_b_n[k] >= 0) 
        model.addConstr(DD_b_n[k] <= M_DD*(1 - α_DD[k])) 
        model.addConstr(DD_b[k] <= M_DD*y_b[k]) 
        
        
    # KCL + DD_b
    for k in Para.L_Node:  # Load Node
        k = int(k)
        model.addConstr(quicksum(Para.A_Matrix[k,b]*(DD_b_p[b] - DD_b_n[b]) for b in range(Para.N_Branch)) == Para.AC_Demand[k] + Para.DC_Demand[k] ) 
    for k in Para.S_Node:  # Substation Node
        k = int(k)    
        model.addConstr(quicksum((Para.A_Matrix[k,b]*(DD_b_p[b] - DD_b_n[b])) for b in range(Para.N_Branch)) == Para.AC_Demand[k] + Para.DC_Demand[k] - Para.GDD_s[k]) 
    
    # DDGC_b
    M_DDGC = Para.GDDGC * 2
    for k in range(Para.N_Branch):
        model.addConstr(DDGC_b[k] == DDGC_b_p[k] + DDGC_b_n[k]) 
        model.addConstr(DDGC_b_p[k] >= 0) 
        model.addConstr(DDGC_b_p[k] <= M_DDGC*α_DDGC[k]) 
        model.addConstr(DDGC_b_n[k] >= 0) 
        model.addConstr(DDGC_b_n[k] <= M_DDGC*(1 - α_DDGC[k])) 
        model.addConstr(DDGC_b[k] <= M_DDGC*y_b[k]) 
    
    # KCL + DDGC_b
    for k in Para.L_Node:  # Load Node
        k = int(k)
        model.addConstr(quicksum(Para.A_Matrix[k,b]*(DDGC_b_p[b] - DDGC_b_n[b]) for b in range(Para.N_Branch)) == Para.AC_DG_Cap[k] + Para.DC_DG_Cap[k] ) 
    
    for k in Para.S_Node:  # Substation Node
        k = int(k)
        model.addConstr(quicksum(Para.A_Matrix[k,b]*(DDGC_b_p[b] - DDGC_b_n[b]) for b in range(Para.N_Branch)) == Para.AC_DG_Cap[k] + Para.DC_DG_Cap[k] - Para.GDDGC_s[k]) 

    # NEENS
    for i in range(Para.N_Node):
        model.addConstr(NEENS[i] == x_type_i[i]*((Para.IC_FaultRate*Para.IC_RepairTime + Para.ACB_FaultRate*Para.ACB_RepairTime + Para.DCB_FaultRate*Para.DCB_RepairTime)*Para.AC_Demand[i] + Para.DCB_FaultRate*Para.DCB_RepairTime*Para.DC_Demand[i] )
                        + (1 - x_type_i[i]) *((Para.RC_FaultRate*Para.RC_RepairTime + Para.ACB_FaultRate*Para.ACB_RepairTime + Para.DCB_FaultRate*Para.DCB_RepairTime)*Para.DC_Demand[i] ) + Para.ACB_FaultRate*Para.ACB_RepairTime*Para.AC_Demand[i] )

    # %%  (2) SAIFI
    model.addConstr(SAIFI == (quicksum(BSAIFI_A[b] for b in range(Para.N_Branch)) 
                          + quicksum(BSAIFI_B[b] for b in range(Para.N_Branch))
                          + quicksum(BSAIFI_C[b] for b in range(Para.N_Branch))
                          + quicksum(NSAIFI[i] for i in range(Para.N_Node))
                          ) / Para.GN)
    # BSAIFI
    M_BSAIFI = Para.GN * 2
    for k in range(Para.N_Branch):
        model.addConstr(BSAIFI_A[k] - Para.ACL_FaultRate[k] * CCN_b[k] <= M_BSAIFI * (1 - ρ_A_b[k]) ) 
        model.addConstr(BSAIFI_A[k] - Para.ACL_FaultRate[k] * CCN_b[k] >= -M_BSAIFI * (1 - ρ_A_b[k]) ) 
        model.addConstr(BSAIFI_A[k]  <= M_BSAIFI * ρ_A_b[k] ) 
        
        model.addConstr(BSAIFI_B[k] - Para.DCL_FaultRate[k] * CCN_b[k] <= M_BSAIFI * (1 - ρ_B_b[k]) ) 
        model.addConstr(BSAIFI_B[k] - Para.DCL_FaultRate[k] * CCN_b[k] >= -M_BSAIFI * (1 - ρ_B_b[k]) ) 
        model.addConstr(BSAIFI_B[k]  <= M_BSAIFI * ρ_B_b[k] ) 
        
        model.addConstr(BSAIFI_C[k] - Para.DCL_FaultRate[k] * CCN_b[k] <= M_BSAIFI * (1 - ρ_C_b[k]) ) 
        model.addConstr(BSAIFI_C[k] - Para.DCL_FaultRate[k] * CCN_b[k] >= -M_BSAIFI * (1 - ρ_C_b[k]) ) 
        model.addConstr(BSAIFI_C[k]  <= M_BSAIFI * ρ_C_b[k] ) 
    
    # CCN
    M_CCN = Para.GN * 2
    for k in range(Para.N_Branch):
        model.addConstr(CCN_b[k] >= 0)
        model.addConstr(CCN_b[k] <= M_CCN * α_CCN[k])
        model.addConstr(CCN_b[k] >= DCN_b[k] - DDGC_b[k]/Para.AP)
        model.addConstr(CCN_b[k] <= DCN_b[k] - DDGC_b[k]/Para.AP + M_CCN * (1- α_CCN[k]))
        
    # DCN
    M_DCN = Para.GN * 2
    for k in range(Para.N_Branch):
        model.addConstr(DCN_b[k] == DCN_b_p[k] + DCN_b_n[k])
        model.addConstr(DCN_b_p[k] >= 0 )
        model.addConstr(DCN_b_p[k] <= M_DCN * α_DCN[k])
        model.addConstr(DCN_b_n[k] >= 0)
        model.addConstr(DCN_b_n[k] <= M_DCN * (1 - α_DCN[k]))
        model.addConstr(DCN_b[k] <= M_DCN * y_b[k])
    
    # KCL + DCN
    for k in Para.L_Node:  # 负荷点
        k = int(k)
        model.addConstr(quicksum(Para.A_Matrix[k,b]*(DCN_b_p[b] - DCN_b_n[b]) for b in range(Para.N_Branch)) == Para.N_Customer[k]) 
    
    for k in Para.S_Node:  # 电源点
        k = int(k)
        model.addConstr(quicksum(Para.A_Matrix[k,b]*(DCN_b_p[b] - DCN_b_n[b]) for b in range(Para.N_Branch)) == Para.N_Customer[k] - Para.GN_s[k]) 
    
    # NSAIFI
    for i in range(Para.N_Node):
        model.addConstr(NSAIFI[i] == x_type_i[i] * ((Para.IC_FaultRate + Para.ACB_FaultRate + Para.DCB_FaultRate) * Para.N_AC_Customer[i] + Para.DCB_FaultRate * Para.N_DC_Customer[i]) + 
                        (1 - x_type_i[i]) * ((Para.RC_FaultRate + Para.ACB_FaultRate + Para.DCB_FaultRate) * Para.N_DC_Customer[i] + Para.ACB_FaultRate * Para.N_AC_Customer[i]))
    
    # %% (3) SAIDI
    model.addConstr(SAIDI == (quicksum(BSAIDI_A[b] for b in range(Para.N_Branch)) 
                          + quicksum(BSAIDI_B[b] for b in range(Para.N_Branch))
                          + quicksum(BSAIDI_C[b] for b in range(Para.N_Branch))
                          + quicksum(NSAIDI[i] for i in range(Para.N_Node))
                          ) / Para.GN)
    # BSAIDI
    M_BSAIDI = 10 * Para.GN
    for k in range(Para.N_Branch):
        model.addConstr(BSAIDI_A[k] - Para.Du_A_b[k] * CCN_b[k] <= M_BSAIDI * (1 - ρ_A_b[k]) ) 
        model.addConstr(BSAIDI_A[k] - Para.Du_A_b[k] * CCN_b[k] >= -M_BSAIDI * (1 - ρ_A_b[k]) ) 
        model.addConstr(BSAIDI_A[k]  <= M_BSAIDI * ρ_A_b[k] ) 
        
        model.addConstr(BSAIDI_B[k] - Para.Du_B_b[k] * CCN_b[k] <= M_BSAIDI * (1 - ρ_B_b[k]) ) 
        model.addConstr(BSAIDI_B[k] - Para.Du_B_b[k] * CCN_b[k] >= -M_BSAIDI * (1 - ρ_B_b[k]) ) 
        model.addConstr(BSAIDI_B[k]  <= M_BSAIDI * ρ_B_b[k] ) 
        
        model.addConstr(BSAIDI_C[k] - Para.Du_C_b[k] * CCN_b[k] <= M_BSAIDI * (1 - ρ_C_b[k]) ) 
        model.addConstr(BSAIDI_C[k] - Para.Du_C_b[k] * CCN_b[k] >= -M_BSAIDI * (1 - ρ_C_b[k]) ) 
        model.addConstr(BSAIDI_C[k]  <= M_BSAIDI * ρ_C_b[k] ) 
    
    # NSAIDI
    for i in range(Para.N_Node):
        model.addConstr(NSAIDI[i] == x_type_i[i] * ((Para.IC_FaultRate * Para.IC_RepairTime + Para.ACB_FaultRate * Para.ACB_RepairTime +  Para.DCB_FaultRate * Para.DCB_RepairTime) * Para.N_AC_Customer[i] + Para.DCB_FaultRate * Para.DCB_RepairTime * Para.N_DC_Customer[i]) + 
                        (1 - x_type_i[i]) * ((Para.RC_FaultRate * Para.RC_RepairTime + Para.ACB_FaultRate * Para.ACB_RepairTime +  Para.DCB_FaultRate * Para.DCB_RepairTime) * Para.N_DC_Customer[i] + Para.ACB_FaultRate * Para.ACB_RepairTime * Para.N_AC_Customer[i]))
    
    # # Reliability Constraint
    # model.addConstr(SAIFI <= 2)
    # model.addConstr(SAIDI <= 8)
    
    # %% Power Flow Equation
    # (1) KVL
    M_V = 6.8*1.05
    for k in range(Para.N_Branch):
        model.addConstr(Para.Z_Branch_A[k] * I_b[k] - quicksum(Para.A_Matrix[i,k] * V_i[i] for i in range(Para.N_Node)) >= -M_V * (1 - ρ_A_b[k]))
        model.addConstr(Para.Z_Branch_A[k] * I_b[k] - quicksum(Para.A_Matrix[i,k] * V_i[i] for i in range(Para.N_Node)) <= M_V * (1 - ρ_A_b[k]))
        model.addConstr(Para.Z_Branch_B[k] * I_b[k] - quicksum((Para.A_Matrix[i,k] * (V_i[i] /Para.VSC_Kc/Para.VSC_Mc + δ_i[i] * (1 - 1/Para.VSC_Kc/Para.VSC_Mc))) for i in range(Para.N_Node)) >= -M_V * (1 - ρ_B_b[k]))
        model.addConstr(Para.Z_Branch_B[k] * I_b[k] - quicksum((Para.A_Matrix[i,k] * (V_i[i] /Para.VSC_Kc/Para.VSC_Mc + δ_i[i] * (1 - 1/Para.VSC_Kc/Para.VSC_Mc))) for i in range(Para.N_Node)) <= M_V * (1 - ρ_B_b[k]))
        model.addConstr(Para.Z_Branch_C[k] * I_b[k] - quicksum(Para.A_Matrix[i,k] * V_i[i] for i in range(Para.N_Node)) >= -M_V * (1 - ρ_C_b[k]))
        model.addConstr(Para.Z_Branch_C[k] * I_b[k] - quicksum(Para.A_Matrix[i,k] * V_i[i] for i in range(Para.N_Node)) <= M_V * (1 - ρ_C_b[k]))
    
    for i in range(Para.N_Node):
        model.addConstr(δ_i[i] >= -M_V * x_type_i[i])
        model.addConstr(δ_i[i] <= M_V * x_type_i[i])
        model.addConstr(δ_i[i] >= - M_V * (1 - x_type_i[i]) + V_i[i])
        model.addConstr(δ_i[i] <= M_V * (1 - x_type_i[i]) + V_i[i])
    
    for i in Para.S_Node:
        i = int(i)
        model.addConstr(V_i[i] == 4.16) # 电源点参考电压
    
    # (2) KCL
    M_I = Para.GDD/10
    for i in Para.L_Node:
        i = int(i)
        model.addConstr(quicksum(Para.A_Matrix[i,b] * I_b[b] for b in range(Para.N_Branch)) +
                         (Para.DC_Demand[i] + Para.AC_Demand[i]/Para.η_INV - Para.AC_DG_Cap[i] 
                          * Para.η_REC - Para.DC_DG_Cap[i]) / Para.Vr_DC >= -M_I * (1-x_type_i[i]))
        model.addConstr(quicksum(Para.A_Matrix[i,b] * I_b[b] for b in range(Para.N_Branch)) +
                         (Para.DC_Demand[i] + Para.AC_Demand[i]/Para.η_INV - Para.AC_DG_Cap[i] 
                          * Para.η_REC - Para.DC_DG_Cap[i]) / Para.Vr_DC <= M_I * (1-x_type_i[i]))
        model.addConstr(quicksum((Para.A_Matrix[i,b] * (I_b[b] + (1 / (3**0.5 * Para.VSC_Kc * Para.VSC_Mc * Para.η_VSC) - 1) * v_b[b])) for b in range(Para.N_Branch)) +
                          (Para.DC_Demand[i]/Para.η_REC + Para.AC_Demand[i] - Para.AC_DG_Cap[i] 
                            - Para.DC_DG_Cap[i] * Para.η_INV) / Para.Vr_AC / 3**0.5 >= -M_I * x_type_i[i])
        model.addConstr(quicksum((Para.A_Matrix[i,b] * (I_b[b] + (1 / (3**0.5 * Para.VSC_Kc * Para.VSC_Mc * Para.η_VSC) - 1) * v_b[b])) for b in range(Para.N_Branch)) +
                          (Para.DC_Demand[i]/Para.η_REC + Para.AC_Demand[i] - Para.AC_DG_Cap[i] 
                            - Para.DC_DG_Cap[i] * Para.η_INV) / Para.Vr_AC / 3**0.5 <= M_I * x_type_i[i])
    
    for i in Para.S_Node:
        i = int(i)
        model.addConstr(quicksum((Para.A_Matrix[i,b] * (I_b[b] + (1 / (3**0.5 * Para.VSC_Kc * Para.VSC_Mc * Para.η_VSC) - 1) * v_b[b])) for b in range(Para.N_Branch)) +
                          (Para.DC_Demand[i]/Para.η_REC + Para.AC_Demand[i] - Para.AC_DG_Cap[i] 
                            - Para.DC_DG_Cap[i] * Para.η_INV) / Para.Vr_AC / 3**0.5 - P_s[i] / Para.Vr_AC / 3**0.5 == 0)

    for k in range(Para.N_Branch):
        model.addConstr(v_b[k] >= - M_I * ρ_B_b[k])
        model.addConstr(v_b[k] <= M_I * ρ_B_b[k])
        model.addConstr(v_b[k] >= - M_I * (1-ρ_B_b[k]) + I_b[k]) 
        model.addConstr(v_b[k] <= M_I * (1-ρ_B_b[k]) + I_b[k])
    
    # %% Constraint: Power Flow Constraint
    # # Nodal Voltage Constraint
    for k in range(Para.N_Node):
        model.addConstr(V_i[k] >= x_type_i[k] * Para.Vmin_DC + (1 - x_type_i[k]) * Para.Vmin_AC)
        model.addConstr(V_i[k] <= x_type_i[k] * Para.Vmax_DC + (1 - x_type_i[k]) * Para.Vmax_AC)
    
    # # Branch Current Constraint
    for k in range(Para.N_Branch):
        model.addConstr(I_b[k] <= ρ_A_b[k] * Para.Imax_AC + (ρ_B_b[k] + ρ_C_b[k]) * Para.Imax_DC)
        model.addConstr(I_b[k] >= -(ρ_A_b[k] * Para.Imax_AC + (ρ_B_b[k] + ρ_C_b[k]) * Para.Imax_DC))
    
    # %% Constraint: Radial Constraint
    # Radial Condition I:
    model.addConstr(quicksum(y_b[b] for b in range(Para.N_Branch)) == Para.N_Node - Para.N_S_Node)
    # # Radial Condition II:
    # for l in range(Para.LoopNum):
    #     model.addConstr(quicksum(y_b[m-1] for m in (Para.Loop[l])) <= len(Para.Loop[l]) - 1)
    

    # %% Fixed Topology (TEST)
    # for k in range(Para.N_Node): # (2)
    #     if k != 0:
    #         model.addConstr(x_type_i[k] == 0) # 
    # for k in range(Para.N_Node): # (3)
    #     if k != 0:
    #         model.addConstr(x_type_i[k] == 1) 
    

    # %% Optimization
    model.setObjective(COST)
    # model.write('model.lp')
    # x_type_i[0].VarHintVal = 0
    # model.setParam('MIPFocus', 3)
    model.setParam('OutputFlag', 1)
    # model.OutputFlag = 0
    model.optimize()
    
    # Debug
    # model.ObjBoundC
    # model.computeIIS()
    # model.write("A.ilp")
    
    # Store result
    Sol = {} 
    x_node = [round(x_type_i[k].x) for k in range(len(x_type_i))]
    Sol['x_node'] = x_node
    y_b = [round(y_b[k].x) for k in range(len(y_b))]
    Sol['y_b_INS'] = y_b
    
    ρ_A_b = [round(ρ_A_b[k].x) for k in range(len(ρ_A_b))]
    ρ_B_b = [round(ρ_B_b[k].x) for k in range(len(ρ_B_b))]
    ρ_C_b = [round(ρ_C_b[k].x) for k in range(len(ρ_C_b))]
    Sol['ρ_B_b'] = ρ_B_b
    Sol['ρ_A_b'] = ρ_A_b
    Sol['ρ_C_b'] = ρ_C_b
    
    V = [V_i[k].x for k in range(Para.N_Node)]
    I = [I_b[b].x * 1000 for b in range(Para.N_Branch)]
    Sol['V'] = V
    Sol['I'] = I
    
    BEENS_A = [BEENS_A[k].x for k in range(Para.N_Branch)]
    BEENS_B = [BEENS_B[k].x for k in range(Para.N_Branch)]
    BEENS_C = [BEENS_C[k].x for k in range(Para.N_Branch)]
    BEENS = np.array(BEENS_A) + np.array(BEENS_B) + np.array(BEENS_C)
    Sol['BEENS'] = BEENS
    
    BSAIDI_A = [BSAIDI_A[k].x for k in range(Para.N_Branch)]
    BSAIDI_B = [BSAIDI_B[k].x for k in range(Para.N_Branch)]
    BSAIDI_C = [BSAIDI_C[k].x for k in range(Para.N_Branch)]
    BSAIDI = np.array(BSAIDI_A) + np.array(BSAIDI_B) + np.array(BSAIDI_C)
    Sol['BSAIDI'] = BSAIDI
    
    BSAIFI_A = [BSAIFI_A[k].x for k in range(Para.N_Branch)]
    BSAIFI_B = [BSAIFI_B[k].x for k in range(Para.N_Branch)]
    BSAIFI_C = [BSAIFI_C[k].x for k in range(Para.N_Branch)]
    BSAIFI = np.array(BSAIFI_A) + np.array(BSAIFI_B) + np.array(BSAIFI_C)
    Sol['BSAIFI'] = BSAIFI
    
    
    NEENS = [NEENS[k].x for k in range(Para.N_Node)]
    NSAIDI = [NSAIDI[k].x for k in range(Para.N_Node)]
    NSAIFI = [NSAIFI[k].x for k in range(Para.N_Node)]
    
    
    # Plot
    PlotFunc.Figure(Para, Sol)
    
    # Output
    ASAI = 1-SAIDI.x/8760
    CAIDI = SAIDI.x/SAIFI.x
    print('\n')
    print('SAIFI：',SAIFI.x)
    print('SAIDI：',SAIDI.x)
    print('EENS：',EENS.x)
    print('ASAI：',ASAI)
    print('CAIDI：',CAIDI)
    
    print('\n')
    print('Investment Cost:',C_I.x)
    print('Maintenance Cost:',C_M.x)
    print('Production Cost:',C_E.x)
    print('Reliability Regulation Cost:',C_Rel.x)
    
    print('\n')
    print('Net Present Cost:',C_Eco.x)
    print('Reliability Cost:',C_Rel.x)
    
    print('\n')
    print('Objective Function：',COST.x)
    
    # %% Save
    # Cost
    result_C = {}
    result_C['Cost_Eco'] = C_Eco.x
    result_C['Cost_Rel'] = C_Rel.x
    result_C['Cost_Inv'] = C_I.x
    result_C['Cost_Mag'] = C_M.x
    result_C['Cost_Ele'] = C_E.x
    result_C = DataFrame(result_C,index=['Cost'])
    
    
    # System Relibility
    result_R = {}
    result_R['SAIFI'] = SAIFI.x
    result_R['SAIDI'] = SAIDI.x
    result_R['EENS'] = EENS.x
    result_R['ASAI'] = ASAI
    result_R['CAIDI'] = CAIDI
    result_R = DataFrame(result_R,index=['Reliability'])
    
    
    # Node Relibility
    result_NRel = {}
    result_NRel['NEENS'] = NEENS
    result_NRel['NSAIDI'] = NSAIDI
    result_NRel['NSAIFI'] = NSAIFI
    result_NRel = DataFrame(result_NRel)
    
    # Branch Reliability
    result_BRel = {}
    
    result_BRel['BEENS_A'] = BEENS_A
    result_BRel['BEENS_B'] = BEENS_B
    result_BRel['BEENS_C'] = BEENS_C
    result_BRel['BEENS'] = BEENS
    
    result_BRel['BSAIDI_A'] = BSAIDI_A
    result_BRel['BSAIDI_B'] = BSAIDI_B
    result_BRel['BSAIDI_C'] = BSAIDI_C
    result_BRel['BSAIDI'] = BSAIDI
    
    result_BRel['BSAIFI_A'] = BSAIFI_A
    result_BRel['BSAIFI_B'] = BSAIFI_B
    result_BRel['BSAIFI_C'] = BSAIFI_C
    result_BRel['BSAIFI'] = BSAIFI
    
    result_BRel = DataFrame(result_BRel)
    
    # Solution
    result_x = {}
    result_y = {}
    result_x['x_node'] = x_node
    result_x['V_node'] = V
    result_y['y_line'] = y_b
    result_y['ρ_A_b'] = ρ_A_b
    result_y['ρ_B_b'] = ρ_B_b
    result_y['ρ_C_b'] = ρ_C_b
    result_y['I_line'] = I
    result_x = DataFrame(result_x)
    result_y = DataFrame(result_y)
    
    # Save to Excel
    # writer = pd.ExcelWriter("(1)Result.xlsx")
    # writer = pd.ExcelWriter("(2)Result.xlsx")
    # writer = pd.ExcelWriter("(3)Result.xlsx")
    # writer = pd.ExcelWriter("(4)Result.xlsx")
    # writer = pd.ExcelWriter("(5)Result.xlsx")
    # writer = pd.ExcelWriter("(6)Result.xlsx")
    
    # writer = pd.ExcelWriter("(1)ComprehensiveOptimization.xlsx")
    # writer = pd.ExcelWriter("(2)ACDN.xlsx")
    # writer = pd.ExcelWriter("(3)DCDN.xlsx")
    # writer = pd.ExcelWriter("(4)Inverse Power Flow Analysis.xlsx")
    
    # result_C.to_excel(writer, sheet_name='result_C')
    # result_R.to_excel(writer, sheet_name='result_R')
    # result_x.to_excel(writer, sheet_name='result_x')
    # result_y.to_excel(writer, sheet_name='result_y')
    # result_BRel.to_excel(writer,sheet_name = 'result_BRel')
    # result_NRel.to_excel(writer,sheet_name = 'result_NRel')
    
    # writer.save()
    