
# coding: utf-8

# In[6]:


import pandas as pd
import numpy as np
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np


# In[7]:


csv_directory = 'C:/Users/SWC/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/results/N2c/20180703_I_N2c_HSVRepeat_EVOS10x/particle_size_500/Triangle/cell_count_GFP.tif.csv'



# In[8]:


df = pd.read_csv(csv_directory)


# In[19]:


part_150 = df[0:16896]
part_250 = df[16896:33792]
part_500 = df[33792:50688]

part_size = ['part_150', 'part_250', 'part_500']
part_size_df = [part_150, part_250, part_500]

# In[20]:

for ab in range(len(part_size)):
    li = part_size_df[ab][0:8448]
    triangle = part_size_df[ab][8448:16896]
    
    li_rfp = li[0:4224]
    li_gfp = li[4224:8448]
    
    triangle_rfp = triangle[0:4224]
    triangle_gfp = triangle[4224:8448]
    
    
    # In[21]:
    
    
    #Flat Condition Means
    Day_list = ["0619", "0621", "0623", "0625", "0627", "0629", "0701", "0703"]
    df_list = [li_rfp, li_gfp, triangle_rfp, triangle_gfp]
    titles = ["li_rfp", "triangle_rfp", "li_gfp", "triangle_gfp"]
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True, sharex=True, figsize = (40,10))
    axes = [ax1, ax2, ax3, ax4]
    
    for aa in range(len(df_list)):
        day_jump = 528
        well_jump = 22
    
    
        US7 = [0]*len(Day_list)
        US8 = [0]*len(Day_list)
        US9 = [0]*len(Day_list)
        US7_8 = [0]*len(Day_list)
        US7_9 = [0]*len(Day_list)
        US8_9 = [0]*len(Day_list)
        US7_8_9 = [0]*len(Day_list)
        GFP = [0]*len(Day_list)
    
        cons = [US7, US8, US9, US7_8, US7_9, US8_9, US7_8_9, GFP]
        legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]
    
        for z in range(len(cons)):
            cond_j = 66*(z+1)
            #print(cond_j)
            for q in range(len(Day_list)):
                cons[z][q] = df_list[aa][q*day_jump:(q+1)*day_jump][cond_j-66:cond_j]["Count"].mean()
                #print(cons[z][q])
    
        for a in range(len(cons)):
            axes[aa].plot(cons[a])
            axes[aa].set_xlabel("Days")
    
        axes[aa].set_title(titles[aa])
        plt.xlabel("Days")
        plt.ylabel("Cell Count")
        axes[aa].legend(legend)
        plt.savefig('C:/Users/SWC/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/HSV_14Jun_Plots/1_Flat_Mean_part_'+part_size[ab] +'.png')
        plt.xticks(np.arange(8))
    
    
    # In[22]:
    
    
    #Flat Condition Mean
    Day_list = ["0619", "0621", "0623", "0625", "0627", "0629", "0701", "0703"]
    df_list = [li_rfp, li_gfp, triangle_rfp, triangle_gfp]
    titles = ["li_rfp", "triangle_rfp", "li_gfp", "triangle_gfp"]
    axes = [ax1, ax2, ax3, ax4]
    
    for aa in range(len(df_list)):
        day_jump = 528
        well_jump = 22
    
        US7 = [0]*len(Day_list)
        US8 = [0]*len(Day_list)
        US9 = [0]*len(Day_list)
        US7_8 = [0]*len(Day_list)
        US7_9 = [0]*len(Day_list)
        US8_9 = [0]*len(Day_list)
        US7_8_9 = [0]*len(Day_list)
        GFP = [0]*len(Day_list)
    
        cons = [US7, US8, US9, US7_8, US7_9, US8_9, US7_8_9, GFP]
        legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]
    
        for z in range(len(cons)):
            cond_j = 66*(z+1)
            #print(cond_j)
            for q in range(len(Day_list)):
                cons[z][q] = df_list[aa][q*day_jump:(q+1)*day_jump][cond_j-66:cond_j]["Count"].mean()
            
            
        fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(1, 8, sharey=True, sharex=True, figsize = (40,10))
        ax1.plot(US7)
        ax1.set_title('US7')
        ax2.plot(US8)
        ax2.set_title('US8')
        ax3.plot(US9)
        ax3.set_title('US9')
        ax4.plot(US7_8)
        ax4.set_title('US7_8')
        ax5.plot(US7_9)
        ax5.set_title('US7_9')
        ax6.plot(US8_9)
        ax6.set_title('US8_9')
        ax7.plot(US7_8_9)
        ax7.set_title('US7_8_9')
        ax8.plot(GFP)
        ax8.set_title('GFP')
        plt.xticks(np.arange(8))
    
    
    # In[23]:
    
    
    #Well mean to condition mean
    
    Day_list = ["0619", "0621", "0623", "0625", "0627", "0629", "0701", "0703"]
    df_list = [li_rfp, li_gfp, triangle_rfp, triangle_gfp]
    titles = ["li_rfp", "li_gfp", "triangle_rfp", "triangle_gfp"]
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True, sharex=True, figsize = (40,10))
    axes = [ax1, ax2, ax3, ax4]
    
    for aa in range(len(df_list)):
        day_jump = 528
        well_jump = 22
    
        US7 = [0]*len(Day_list)
        US8 = [0]*len(Day_list)
        US9 = [0]*len(Day_list)
        US7_8 = [0]*len(Day_list)
        US7_9 = [0]*len(Day_list)
        US8_9 = [0]*len(Day_list)
        US7_8_9 = [0]*len(Day_list)
        GFP = [0]*len(Day_list)
    
        cons = [US7, US8, US9, US7_8, US7_9, US8_9, US7_8_9, GFP]
        legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]
    
    
    
        for v in range(len(cons)):
            well_j = 22*(v+1)
            for q in range(len(Day_list)):
                w1 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j-22:well_j]["Count"].mean()
                w2 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j:well_j+22]["Count"].mean()
                w3 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j+22:well_j+44]["Count"].mean()
                w_m = [w1, w2, w3]
                cons[v][q] = np.mean(w_m)
    
        for a in range(len(cons)):
            axes[aa].plot(cons[a])
            axes[aa].set_xlabel("Days")
        axes[aa].set_title(titles[aa])
        plt.xlabel("Days")
        plt.ylabel("Cell Count")
        axes[aa].legend(legend)
        plt.savefig('C:/Users/SWC/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/HSV_14Jun_Plots/2_Well_Mean_'+part_size[ab] +'.png')
        plt.xticks(np.arange(8))
    
    
    # In[9]:
    
    
    #Well mean to condition mean
    Day_list = ["0619", "0621", "0623", "0625", "0627", "0629", "0701", "0703"]
    df_list = [li_rfp, li_gfp, triangle_rfp, triangle_gfp]
    titles = ["li_rfp", "triangle_rfp", "li_gfp", "triangle_gfp"]
    axes = [ax1, ax2, ax3, ax4]
    
    for aa in range(len(df_list)):
        day_jump = 528
        well_jump = 22
    
        US7 = [0]*len(Day_list)
        US8 = [0]*len(Day_list)
        US9 = [0]*len(Day_list)
        US7_8 = [0]*len(Day_list)
        US7_9 = [0]*len(Day_list)
        US8_9 = [0]*len(Day_list)
        US7_8_9 = [0]*len(Day_list)
        GFP = [0]*len(Day_list)
    
        cons = [US7, US8, US9, US7_8, US7_9, US8_9, US7_8_9, GFP]
        legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]
    
        for v in range(len(cons)):
            well_j = 22*(v+1)
            for q in range(len(Day_list)):
                w1 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j-22:well_j]["Count"].mean()
                w2 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j:well_j+22]["Count"].mean()
                w3 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j+22:well_j+44]["Count"].mean()
                w_m = [w1, w2, w3]
                cons[v][q] = np.mean(w_m)
            
        fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(1, 8, sharey=True, sharex=True, figsize = (40,10))
        ax1.plot(US7)
        ax1.set_title('US7')
        ax2.plot(US8)
        ax2.set_title('US8')
        ax3.plot(US9)
        ax3.set_title('US9')
        ax4.plot(US7_8)
        ax4.set_title('US7_8')
        ax5.plot(US7_9)
        ax5.set_title('US7_9')
        ax6.plot(US8_9)
        ax6.set_title('US8_9')
        ax7.plot(US7_8_9)
        ax7.set_title('US7_8_9')
        ax8.plot(GFP)
        ax8.set_title('GFP')
        plt.xticks(np.arange(8))
    
    
    # In[25]:
    
    
    #Combined Means
    triangle_rfp_2 = triangle_rfp.reset_index(drop=True)
    com_rfp = pd.concat((li_rfp, triangle_rfp_2))
    by_row_rfp = com_rfp.groupby(com_rfp.index)
    com_rfp_means = by_row_rfp.mean()
    
    triangle_gfp_2 = triangle_gfp.reset_index(drop=True)
    com_gfp = pd.concat((li_gfp, triangle_gfp_2))
    by_row_gfp = com_gfp.groupby(com_gfp.index)
    com_gfp_means = by_row_gfp.mean()
    
    
    Day_list = ["0619", "0621", "0623", "0625", "0627", "0629", "0701", "0703"]
    df_list = [com_rfp_means, com_gfp]
    titles = ["com_rfp", "com_gfp"]
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, sharex=True, figsize = (40,10))
    axes = [ax1, ax2]
    
    for aa in range(len(df_list)):
        day_jump = 528
        well_jump = 22
    
    
        US7 = [0]*len(Day_list)
        US8 = [0]*len(Day_list)
        US9 = [0]*len(Day_list)
        US7_8 = [0]*len(Day_list)
        US7_9 = [0]*len(Day_list)
        US8_9 = [0]*len(Day_list)
        US7_8_9 = [0]*len(Day_list)
        GFP = [0]*len(Day_list)
    
        cons = [US7, US8, US9, US7_8, US7_9, US8_9, US7_8_9, GFP]
        legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]
    
        for z in range(len(cons)):
            cond_j = 66*(z+1)
            #print(cond_j)
            for q in range(len(Day_list)):
                cons[z][q] = df_list[aa][q*day_jump:(q+1)*day_jump][cond_j-66:cond_j]["Count"].mean()
                #print(cons[z][q])
    
        for a in range(len(cons)):
            axes[aa].plot(cons[a])
            axes[aa].set_xlabel("Days")
    
        axes[aa].set_title(titles[aa])
        plt.xlabel("Days")
        plt.ylabel("Cell Count")
        axes[aa].legend(legend)
        plt.savefig('C:/Users/SWC/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/HSV_14Jun_Plots/3_Combined_Mean_'+part_size[ab] +'.png')
        plt.xticks(np.arange(8))
    
    
    # In[26]:
    
    
    #Combined Means Percentage Increase
    
    Day_list = ["0619", "0621", "0623", "0625", "0627", "0629", "0701", "0703"]
    df_list = [com_rfp_means, com_gfp_means]
    titles = ["com_rfp", "com_gfp"]
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True,  figsize = (40,10))
    axes = [ax1, ax2, ax3, ax4]
    
    for aa in range(len(df_list)):
        day_jump = 528
        well_jump = 22
    
    
        US7 = [0]*len(Day_list)
        US8 = [0]*len(Day_list)
        US9 = [0]*len(Day_list)
        US7_8 = [0]*len(Day_list)
        US7_9 = [0]*len(Day_list)
        US8_9 = [0]*len(Day_list)
        US7_8_9 = [0]*len(Day_list)
        GFP = [0]*len(Day_list)
        
        perc_US7 = [0]*len(Day_list)
        perc_US8 = [0]*len(Day_list)
        perc_US9 = [0]*len(Day_list)
        perc_US7_8 = [0]*len(Day_list)
        perc_US7_9 = [0]*len(Day_list)
        perc_US8_9 = [0]*len(Day_list)
        perc_US7_8_9 = [0]*len(Day_list)
        perc_GFP = [0]*len(Day_list)
        
        std = [0]*len(Day_list)
    
        cons = [US7, US8, US9, US7_8, US7_9, US8_9, US7_8_9, GFP]
        perc_cons = [perc_US7, perc_US8, perc_US9, perc_US7_8, perc_US7_9, perc_US8_9, perc_US7_8_9, perc_GFP]
        legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]
    
        for z in range(len(cons)):
            cond_j = 66*(z+1)
            #print(cond_j)
            for q in range(len(Day_list)):
                cons[z][q] = df_list[aa][q*day_jump:(q+1)*day_jump][cond_j-66:cond_j]["Count"].mean()
                perc_cons[z][q] = ((cons[z][q]-cons[z][0])/cons[z][0])*100
                #print(perc_cons)
    
        for a in range(len(cons)):
            axes[aa].plot(perc_cons[a])
            axes[aa].set_xlabel("Days")
    
        axes[aa].set_title(titles[aa])
        plt.xlabel("Days")
        plt.ylabel("Cell Count")
        axes[aa].legend(legend)
        plt.xticks(np.arange(8))
        
        
        legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]
        d1m = [perc_US7[4],perc_US8[4],perc_US9[4],perc_US7_8[4],perc_US7_9[4],perc_US8_9[4],perc_US7_8_9[4],perc_GFP[4]]
    
        axes[aa+2].bar(np.arange(len(legend)), d1m, color=['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'grey'])
        axes[aa+2].set_xticks(np.arange(len(legend)))
        axes[aa+2].set_xticklabels(legend, rotation =60)
        axes[aa+2].set_title(titles[aa])
        axes[aa+2].set_xlabel("Condition")
    
        
        
    plt.savefig('C:/Users/SWC/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/HSV_14Jun_Plots/4_Combined_Mean_Perc_Inc_'+part_size[ab] +'.png')
    
    
    # In[27]:
    
    
    #Well Mean Percentage Increase
    
    Day_list = ["0619", "0621", "0623", "0625", "0627", "0629", "0701", "0703"]
    df_list = [com_rfp_means, com_gfp_means]
    titles = ["com_rfp", "com_gfp"]
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True,  figsize = (40,10))
    axes = [ax1, ax2, ax3, ax4]
    
    for aa in range(len(df_list)):
        day_jump = 528
        well_jump = 22
    
    
        US7 = [0]*len(Day_list)
        US8 = [0]*len(Day_list)
        US9 = [0]*len(Day_list)
        US7_8 = [0]*len(Day_list)
        US7_9 = [0]*len(Day_list)
        US8_9 = [0]*len(Day_list)
        US7_8_9 = [0]*len(Day_list)
        GFP = [0]*len(Day_list)
        
        perc_US7 = [0]*len(Day_list)
        perc_US8 = [0]*len(Day_list)
        perc_US9 = [0]*len(Day_list)
        perc_US7_8 = [0]*len(Day_list)
        perc_US7_9 = [0]*len(Day_list)
        perc_US8_9 = [0]*len(Day_list)
        perc_US7_8_9 = [0]*len(Day_list)
        perc_GFP = [0]*len(Day_list)
        
        std = [0]*len(Day_list)
    
        cons = [US7, US8, US9, US7_8, US7_9, US8_9, US7_8_9, GFP]
        perc_cons = [perc_US7, perc_US8, perc_US9, perc_US7_8, perc_US7_9, perc_US8_9, perc_US7_8_9, perc_GFP]
        legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]
    
        for v in range(len(cons)):
            well_j = 22*(v+1)
            for q in range(len(Day_list)):
                w1 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j-22:well_j]["Count"].mean()
                w2 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j:well_j+22]["Count"].mean()
                w3 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j+22:well_j+44]["Count"].mean()
                w_m = [w1, w2, w3]
                cons[v][q] = np.mean(w_m)
                perc_cons[v][q] = ((cons[v][q]-cons[v][0])/cons[v][0])*100
    
        for a in range(len(cons)):
            axes[aa].plot(perc_cons[a])
            axes[aa].set_xlabel("Days")
    
        axes[aa].set_title(titles[aa])
        plt.xlabel("Days")
        plt.ylabel("Cell Count")
        axes[aa].legend(legend)
        plt.xticks(np.arange(8))
        
        
        legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]
        d1m = [perc_US7[4],perc_US8[4],perc_US9[4],perc_US7_8[4],perc_US7_9[4],perc_US8_9[4],perc_US7_8_9[4],perc_GFP[4]]
    
        axes[aa+2].bar(np.arange(len(legend)), d1m, color=['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'grey'])
        axes[aa+2].set_xticks(np.arange(len(legend)))
        axes[aa+2].set_xticklabels(legend, rotation =60)
        axes[aa+2].set_title(titles[aa])
        axes[aa+2].set_xlabel("Condition")
    
        
        
    plt.savefig('C:/Users/SWC/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/HSV_14Jun_Plots/7_Well_Mean_Perc_Inc_'+part_size[ab] +'.png')
    
    
    # In[16]:
    
    
    #Li Means Percentage Increase
    
    Day_list = ["0619", "0621", "0623", "0625", "0627", "0629", "0701", "0703"]
    df_list = [li_rfp, li_gfp]
    titles = ["li_rfp", "li_gfp"]
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True, figsize = (40,10))
    axes = [ax1, ax2, ax3, ax4]
    
    for aa in range(len(df_list)):
        day_jump = 528
        well_jump = 22
    
    
        US7 = [0]*len(Day_list)
        US8 = [0]*len(Day_list)
        US9 = [0]*len(Day_list)
        US7_8 = [0]*len(Day_list)
        US7_9 = [0]*len(Day_list)
        US8_9 = [0]*len(Day_list)
        US7_8_9 = [0]*len(Day_list)
        GFP = [0]*len(Day_list)
        
        perc_US7 = [0]*len(Day_list)
        perc_US8 = [0]*len(Day_list)
        perc_US9 = [0]*len(Day_list)
        perc_US7_8 = [0]*len(Day_list)
        perc_US7_9 = [0]*len(Day_list)
        perc_US8_9 = [0]*len(Day_list)
        perc_US7_8_9 = [0]*len(Day_list)
        perc_GFP = [0]*len(Day_list)
        
        std = [0]*len(Day_list)
    
        cons = [US7, US8, US9, US7_8, US7_9, US8_9, US7_8_9, GFP]
        perc_cons = [perc_US7, perc_US8, perc_US9, perc_US7_8, perc_US7_9, perc_US8_9, perc_US7_8_9, perc_GFP]
        legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]
    
        for v in range(len(cons)):
            well_j = 22*(v+1)
            for q in range(len(Day_list)):
                w1 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j-22:well_j]["Count"].mean()
                w2 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j:well_j+22]["Count"].mean()
                w3 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j+22:well_j+44]["Count"].mean()
                w_m = [w1, w2, w3]
                cons[v][q] = np.mean(w_m)
                perc_cons[v][q] = ((cons[v][q]-cons[v][0])/cons[v][0])*100
    
        for a in range(len(cons)):
            axes[aa].plot(perc_cons[a])
            axes[aa].set_xlabel("Days")
        axes[aa].set_title(titles[aa])
        plt.xlabel("Days")
        plt.ylabel("Cell Count")
        axes[aa].legend(legend)
        plt.xticks(np.arange(8))
        
        
        legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]
        d1m = [perc_US7[4],perc_US8[4],perc_US9[4],perc_US7_8[4],perc_US7_9[4],perc_US8_9[4],perc_US7_8_9[4],perc_GFP[4]]
        
           
        axes[aa+2].bar(np.arange(len(legend)), d1m, color=['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'grey'])
        axes[aa+2].set_xticks(np.arange(len(legend)))
        axes[aa+2].set_xticklabels(legend, rotation =60)
        axes[aa+2].set_title(titles[aa])
        axes[aa+2].set_xlabel("Condition")
    
        
        
    plt.savefig('C:/Users/SWC/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/HSV_14Jun_Plots/5_Li_Mean_Perc_Inc'+part_size[ab] +'.png')
    
    
    # In[28]:
    
    
    #Li and Triangle Means Percentage Increase
    
    Day_list = ["0619", "0621", "0623", "0625", "0627", "0629", "0701", "0703"]
    df_list = [li_rfp, li_gfp, triangle_rfp, triangle_gfp]
    titles = ["li_rfp", "li_gfp", "triangle_rfp", "triangle_gfp"]
    fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(1, 8, sharey=True, figsize = (40,10))
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
    
    for aa in range(len(df_list)):
        day_jump = 528
        well_jump = 22
    
    
        US7 = [0]*len(Day_list)
        US8 = [0]*len(Day_list)
        US9 = [0]*len(Day_list)
        US7_8 = [0]*len(Day_list)
        US7_9 = [0]*len(Day_list)
        US8_9 = [0]*len(Day_list)
        US7_8_9 = [0]*len(Day_list)
        GFP = [0]*len(Day_list)
        
        perc_US7 = [0]*len(Day_list)
        perc_US8 = [0]*len(Day_list)
        perc_US9 = [0]*len(Day_list)
        perc_US7_8 = [0]*len(Day_list)
        perc_US7_9 = [0]*len(Day_list)
        perc_US8_9 = [0]*len(Day_list)
        perc_US7_8_9 = [0]*len(Day_list)
        perc_GFP = [0]*len(Day_list)
        
        std = [0]*len(Day_list)
    
        cons = [US7, US8, US9, US7_8, US7_9, US8_9, US7_8_9, GFP]
        perc_cons = [perc_US7, perc_US8, perc_US9, perc_US7_8, perc_US7_9, perc_US8_9, perc_US7_8_9, perc_GFP]
        legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]
    
    
        for v in range(len(cons)):
            well_j = 22*(v+1)
            for q in range(len(Day_list)):
                w1 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j-22:well_j]["Count"].mean()
                w2 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j:well_j+22]["Count"].mean()
                w3 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j+22:well_j+44]["Count"].mean()
                w_m = [w1, w2, w3]
                cons[v][q] = np.mean(w_m)
                perc_cons[v][q] = ((cons[v][q]-cons[v][0])/cons[v][0])*100
    
    
        for a in range(len(cons)):
            axes[aa].plot(perc_cons[a])
            axes[aa].set_xticks(np.arange(8))
            axes[aa].set_xticklabels(np.arange(8))
            axes[aa].set_xlabel("Days")
        axes[aa].set_title(titles[aa])
        plt.ylabel("% Increase")
        axes[aa].legend(legend)
        
        
        legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]
        d1m = [perc_US7[4],perc_US8[4],perc_US9[4],perc_US7_8[4],perc_US7_9[4],perc_US8_9[4],perc_US7_8_9[4],perc_GFP[4]]
    
        
        
        axes[aa+4].bar(np.arange(len(legend)), d1m, color=['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'grey'])
        axes[aa+4].set_xticks(np.arange(len(legend)))
        axes[aa+4].set_xticklabels(legend, rotation =60)
        axes[aa+4].set_title(titles[aa])
        #axes[aa+4].set_facecolor()
        
    
    
        
    plt.savefig('C:/Users/SWC/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/HSV_14Jun_Plots/6_Li_Triange_Mean_Perc_Inc_part_150'+part_size[ab] +'.png')
    
