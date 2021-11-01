import matplotlib.pyplot as plt
import numpy as np

def Figure(Para,Sol):
    x = Para.Node_x
    y = Para.Node_y
    plt.figure()
    
    # %% Power Flow Figure
    # 线路绘制 
    for n in range(Para.N_Branch):  # 遍历支路
        if Sol['y_b_INS'][n] == 1: # 规划建设的线路
            if Sol['ρ_A_b'][n] == 1: # 交流线路
                x1 = x[int(round(Para.Branch[n,1]))]
                y1 = y[int(round(Para.Branch[n,1]))]
                x2 = x[int(round(Para.Branch[n,2]))]
                y2 = y[int(round(Para.Branch[n,2]))]
                plt.plot([x1,x2],[y1,y2],'r-',linewidth=2)
                Current = round(Sol['I'][n],2)
                plt.text((x1+x2)/2, (y1+y2)/2, '%s'%Current)
            if (Sol['ρ_B_b'][n] == 1) or (Sol['ρ_C_b'][n] == 1)  : # 直流线路
                x1 = x[int(round(Para.Branch[n,1]))]
                y1 = y[int(round(Para.Branch[n,1]))]
                x2 = x[int(round(Para.Branch[n,2]))]
                y2 = y[int(round(Para.Branch[n,2]))]
                plt.plot([x1,x2],[y1,y2],'b-',linewidth=2)
                Current = round(Sol['I'][n],2)
                plt.text((x1+x2)/2, (y1+y2)/2, '%s'%Current)
        else: # 不建设的线路
            x1 = x[int(round(Para.Branch[n,1]))]
            y1 = y[int(round(Para.Branch[n,1]))]
            x2 = x[int(round(Para.Branch[n,2]))]
            y2 = y[int(round(Para.Branch[n,2]))]
            plt.plot([x1,x2],[y1,y2],'k--',linewidth=0.5)
            Current = round(Sol['I'][n],2)
            plt.text((x1+x2)/2, (y1+y2)/2, '%s'%Current)
    # 节点绘制
    for n in range(Para.N_Node):  
        if Sol['x_node'][n] == 1: # 直流节点
            plt.plot(x[n],y[n],'bs',markersize=10)
        else: # 交流节点
            plt.plot(x[n],y[n],'rs',markersize=10)
            
    # 节点绘制  Bus
    for n in range(Para.N_Node):  
        # plt.text(x[n], y[n] + 0.03, '%s'%n)
        Vol = round(Sol['V'][n],2)
        plt.text(x[n], y[n] + 0.03, '%s'%Vol)
        
        # plt.plot(x[n],y[n],'rs',markersize=10)
        # plt.plot(x[n],y[n],'rs',markersize=10)

        
    plt.axis('equal')
    # plt.title('Power Flow Result')
    plt.show()
    
    # # %% Reliability Figure
    # plt.figure()
    # # 线路绘制 
    # for n in range(Para.N_Branch):  # 遍历支路
    #     if Sol['y_b_INS'][n] == 1: # 规划建设的线路
    #         if Sol['ρ_A_b'][n] == 1: # 交流线路
    #             x1 = x[int(round(Para.Branch[n,1]))]
    #             y1 = y[int(round(Para.Branch[n,1]))]
    #             x2 = x[int(round(Para.Branch[n,2]))]
    #             y2 = y[int(round(Para.Branch[n,2]))]
    #             plt.plot([x1,x2],[y1,y2],'r-',linewidth=2)
    #             Current = round(Sol['I'][n],2) # BEENS
    #             plt.text((x1+x2)/2, (y1+y2)/2, '%s'%Current)
    #         if (Sol['ρ_B_b'][n] == 1) or (Sol['ρ_C_b'][n] == 1)  : # 直流线路
    #             x1 = x[int(round(Para.Branch[n,1]))]
    #             y1 = y[int(round(Para.Branch[n,1]))]
    #             x2 = x[int(round(Para.Branch[n,2]))]
    #             y2 = y[int(round(Para.Branch[n,2]))]
    #             plt.plot([x1,x2],[y1,y2],'b-',linewidth=2)
    #             Current = round(Sol['I'][n],2) # BEENS
    #             plt.text((x1+x2)/2, (y1+y2)/2, '%s'%Current)
    #     else: # 不建设的线路
    #         x1 = x[int(round(Para.Branch[n,1]))]
    #         y1 = y[int(round(Para.Branch[n,1]))]
    #         x2 = x[int(round(Para.Branch[n,2]))]
    #         y2 = y[int(round(Para.Branch[n,2]))]
    #         plt.plot([x1,x2],[y1,y2],'k--',linewidth=0.5)
    #         Current = round(Sol['I'][n],2) # BEENS
    #         plt.text((x1+x2)/2, (y1+y2)/2, '%s'%Current)
    # # 节点绘制
    # for n in range(Para.N_Node):  
    #     if Sol['x_node'][n] == 1: # 直流节点
    #         plt.plot(x[n],y[n],'bs',markersize=10)
    #     else: # 交流节点
    #         plt.plot(x[n],y[n],'rs',markersize=10)
            
    # # 节点绘制  Bus
    # for n in range(Para.N_Node):  
    #     # plt.text(x[n], y[n] + 0.03, '%s'%n)
    #     Vol = round(Sol['V'][n],2)
    #     plt.text(x[n], y[n] + 0.03, '%s'%Vol)
        
    #     # plt.plot(x[n],y[n],'rs',markersize=10)
    #     # plt.plot(x[n],y[n],'rs',markersize=10)

    
    # plt.axis('equal')
    # plt.title('Reliability Figure')
    # plt.show()