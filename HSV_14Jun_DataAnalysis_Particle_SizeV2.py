
# coding: utf-8

# In[246]:


import pandas as pd
import numpy as np
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np

import numpy as np
import scipy.stats

csv_directory = 'C:/Users/SimonThompson/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/results/20180703_I_N2c_HSVRepeat_EVOS10x/particle_size_3000/Triangle/cell_count_GFP.tif.csv'

df = pd.read_csv(csv_directory)


part_1500 = df[0:16896]
part_3000 = df[16896:33792]

part_size = ['part_1500', 'part_3000']
part_size_df = [part_1500, part_3000]





for ab in range(len(part_size)):
    li = part_size_df[ab][0:8448]
    triangle = part_size_df[ab][8448:16896]
    
    li_rfp = li[0:4224]
    li_gfp = li[4224:8448]
    
    triangle_rfp = triangle[0:4224]
    triangle_gfp = triangle[4224:8448]
    
    triangle_rfp_2 = triangle_rfp.reset_index(drop=True)
    com_rfp = pd.concat((li_rfp, triangle_rfp_2))
    by_row_rfp = com_rfp.groupby(com_rfp.index)
    com_rfp_means = by_row_rfp.mean()

    triangle_gfp_2 = triangle_gfp.reset_index(drop=True)
    com_gfp = pd.concat((li_gfp, triangle_gfp_2))
    by_row_gfp = com_gfp.groupby(com_gfp.index)
    com_gfp_means = by_row_gfp.mean()




    #  1  Mean across whole condition, 66 images, total cell count

    Day_list = ["0619", "0621", "0623", "0625", "0627", "0629", "0701", "0703"]
    df_list = [li_rfp, li_gfp, triangle_rfp, triangle_gfp]
    titles = ["li_rfp", "triangle_rfp", "li_gfp", "triangle_gfp"]
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True, sharex=True, figsize = (40,10))
    axes = [ax1, ax2, ax3, ax4]
    legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]

    for aa in range(len(df_list)):
        day_jump = 528
        well_jump = 22

        for z in range(len(legend)):
            cond_j = 66*(z+1)
            c1 = [0]*len(Day_list)
            c1_sem = [0]*len(Day_list)

            for q in range(len(Day_list)):
                c1[q] = df_list[aa][q*day_jump:(q+1)*day_jump][cond_j-66:cond_j]["Count"].mean()
                c1_sem[q] = df_list[aa][q*day_jump:(q+1)*day_jump][cond_j-66:cond_j]["Count"].std()

            axes[aa].errorbar(Day_list, c1, yerr = c1_sem)

        axes[aa].set_xlabel("Days")
        axes[aa].set_ylabel("Cell Count Per Image")
        axes[aa].set_title(titles[aa])
        axes[aa].legend(legend)    
        plt.savefig('C:/Users/SimonThompson/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/HSV_14Jun_Plots/1_Flat_Mean_Total_Cell_'+ part_size[ab] + '.png')


    # In[291]:


    #2 Means across wells and then averaged, total cell count
    Day_list = ["0619", "0621", "0623", "0625", "0627", "0629", "0701", "0703"]
    df_list = [li_rfp, li_gfp, triangle_rfp, triangle_gfp]
    titles = ["li_rfp", "li_gfp", "triangle_rfp", "triangle_gfp"]
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True, sharex=True, figsize = (40,10))
    axes = [ax1, ax2, ax3, ax4]

    legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]

    for aa in range(len(df_list)):
        day_jump = 528
        well_jump = 66

        for v in range(len(legend)):        
            well_j = 66*(v+1)
            c1 = [0]*len(Day_list)
            c1_sem = [0]*len(Day_list)

            for q in range(len(Day_list)):
                w1 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j-well_jump:well_j-int(((well_jump/3)*2))]["Count"].mean()
                w2 = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)*2):int(well_j-(well_jump/3))]["Count"].mean()
                w3 = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)):well_j]["Count"].mean()
                w_m = [w1, w2, w3]
                c1[q] = np.mean(w_m)

                w1_sem = df_list[aa][q*day_jump:(q+1)*day_jump][well_j-well_jump:well_j-int(((well_jump/3)*2))]["Count"].sem()
                w2_sem = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)*2):int(well_j-(well_jump/3))]["Count"].sem()
                w3_sem = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)):well_j]["Count"].sem()



            axes[aa].plot(Day_list, c1)
            axes[aa].set_xlabel("Days")

        axes[aa].set_title(titles[aa])
        axes[aa].set_ylabel("Cell Count")
        axes[aa].legend(legend)

    plt.savefig('C:/Users/SimonThompson/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/HSV_14Jun_Plots/2_Well_Means_Total_Cell_' + part_size[ab]+  '.png')


    # In[353]:


    #3 Combined Means of Li and Triangle across whole condition, total cell count
    triangle_rfp_2 = triangle_rfp.reset_index(drop=True)
    com_rfp = pd.concat((li_rfp, triangle_rfp_2))
    by_row_rfp = com_rfp.groupby(com_rfp.index)
    com_rfp_means = by_row_rfp.mean()

    Day_list = ["0619", "0621", "0623", "0625", "0627", "0629", "0701", "0703"]
    df_list = [com_rfp_means, com_gfp_means]
    titles = ["com_rfp", "com_gfp"]
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, sharex=True, figsize = (40,10))
    axes = [ax1, ax2]
    legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]

    for aa in range(len(df_list)):
        day_jump = 528
        for z in range(len(legend)):
            cond_j = 66*(z+1)
            c1 = [0]*len(Day_list)
            c1_sem = [0]*len(Day_list)


            for q in range(len(Day_list)):
                c1[q] = df_list[aa][q*day_jump:(q+1)*day_jump][cond_j-66:cond_j]["Count"].mean()
                c1_sem[q] = df_list[aa][q*day_jump:(q+1)*day_jump][cond_j-66:cond_j]["Count"].std()


            axes[aa].errorbar(Day_list, c1, yerr = c1_sem)
            axes[aa].set_xlabel("Days")
            axes[aa].set_ylabel("Cell Count")

        axes[aa].set_title(titles[aa])
        axes[aa].legend(legend)

    plt.savefig('C:/Users/SimonThompson/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/HSV_14Jun_Plots/3_Combined_Means_Li_Triangle_Total_Cell_'+ part_size[ab] + '.png')


    # In[293]:


    #4 Combined Means of Li and Triangle across wells, total cell count
    Day_list = ["0619", "0621", "0623", "0625", "0627", "0629", "0701", "0703"]
    df_list = [com_rfp_means, com_gfp_means]
    titles = ["com_rfp", "com_gfp"]
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,  figsize = (40,10))
    axes = [ax1, ax2]
    legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]

    for aa in range(len(df_list)):
        day_jump = 528
        well_jump = 66

        for v in range(len(legend)):
            well_j = 66*(v+1)
            c1 = [0]*len(Day_list)

            for q in range(len(Day_list)):
                w1 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j-well_jump:well_j-int(((well_jump/3)*2))]["Count"].mean()
                w2 = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)*2):int(well_j-(well_jump/3))]["Count"].mean()
                w3 = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)):well_j]["Count"].mean()
                w_m = [w1, w2, w3]
                c1[q] = np.mean(w_m)

            axes[aa].plot(Day_list, c1)
            axes[aa].set_xlabel("Days")
            axes[aa].set_ylabel("Cell Count")

        axes[aa].set_title(titles[aa])
        axes[aa].legend(legend)

    plt.savefig('C:/Users/SimonThompson/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/HSV_14Jun_Plots/4_Combined_Mean_Li_Triangle_Well_Total_Cell_'+ part_size[ab] + '.png')


    # In[372]:


    #5 Combined (Li + Triangle), Well mean, well comparison, total cell count
    triangle_gfp_2 = triangle_gfp.reset_index(drop=True)
    com_gfp = pd.concat((li_gfp, triangle_gfp_2))
    by_row_gfp = com_gfp.groupby(com_gfp.index)
    com_gfp_means = by_row_gfp.mean()

    Day_list = ["0619", "0621", "0623", "0625", "0627", "0629", "0701", "0703"]
    df_list = [com_rfp_means]
    titles = ["com_rfp"]
    fig, ((ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8)) = plt.subplots(1, 8, sharey=True, sharex=True, figsize = (40,10))
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

    for aa in range(len(df_list)):
        day_jump = 528
        well_jump = 66

        legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]

        for v in range(len(legend)):
            well_j = 66*(v+1)

            w1 = [0]*len(Day_list)
            w2 = [0]*len(Day_list)
            w3 = [0]*len(Day_list)
            w_con = [0]*len(Day_list)

            w1_sem = [0]*len(Day_list)
            w2_sem = [0]*len(Day_list)
            w3_sem = [0]*len(Day_list)
            w_con_sem = [0]*len(Day_list)

            for q in range(len(Day_list)): 
                w1[q] = df_list[aa][q*day_jump:(q+1)*day_jump][well_j-well_jump:well_j-int(((well_jump/3)*2))]["Count"].mean()
                w2[q] = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)*2):int(well_j-(well_jump/3))]["Count"].mean()
                w3[q] = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)):well_j]["Count"].mean()

                w_con[q] = np.mean([w1[q], w2[q], w3[q]])

                w1_sem[q] = df_list[aa][q*day_jump:(q+1)*day_jump][well_j-well_jump:well_j-int(((well_jump/3)*2))]["Count"].std()
                w2_sem[q] = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)*2):int(well_j-(well_jump/3))]["Count"].std()
                w3_sem[q] = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)):well_j]["Count"].std()


            axes[v + (aa*8)].errorbar(Day_list, w1, label = 'W1', yerr = w1_sem)
            axes[v + (aa*8)].errorbar(Day_list, w2, label = 'W2', yerr = w2_sem)
            axes[v + (aa*8)].errorbar(Day_list, w3, label = 'W3', yerr = w3_sem)
            axes[v + (aa*8)].plot(Day_list, w_con, label = 'Mean')
            axes[v + (aa*8)].legend()
            axes[v + (aa*8)].set_xlabel("Days")
            axes[v + (aa*8)].set_title(legend[v])

    plt.savefig('C:/Users/SimonThompson/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/HSV_14Jun_Plots/5_Well_Mean_Well_Comp_Total_Cell_'+ part_size[ab] + '.png')


    # In[373]:


    #6 Combined Means (Li + Triangle),whole condition, percentage Increase

    Day_list = ["0619", "0621", "0623", "0625", "0627", "0629", "0701", "0703"]
    df_list = [com_rfp_means, com_gfp_means]
    titles = ["com_rfp", "com_gfp"]
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True,  figsize = (40,10))
    axes = [ax1, ax2, ax3, ax4]
    legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]

    for aa in range(len(df_list)):
        day_jump = 528
        well_jump = 22   
        bar = [0]*len(legend)
        bar_sem = [0]*len(legend)

        for z in range(len(legend)):
            cond_j = 66*(z+1)
            c1 = [0]*len(Day_list)
            perc_c1 = [0]*len(Day_list)
            c1_sem = [0]*len(Day_list)

            for q in range(len(Day_list)):
                c1[q] = df_list[aa][q*day_jump:(q+1)*day_jump][cond_j-66:cond_j]["Count"].mean()
                perc_c1[q] = ((c1[q]-c1[0])/c1[0])*100
                c1_sem[q] = df_list[aa][q*day_jump:(q+1)*day_jump][cond_j-66:cond_j]["Count"].sem()

                if q == 4:
                    bar[z] = perc_c1[q]
                    bar_sem[z] = c1_sem[q]

            axes[aa].errorbar(Day_list, perc_c1)
            axes[aa].set_xlabel("Days")
            axes[aa].set_ylabel("% Increase")

        axes[aa+2].bar(np.arange(len(legend)), bar ,  color=['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'grey'])
        axes[aa+2].set_xticks(np.arange(len(legend)))
        axes[aa+2].set_xticklabels(legend, rotation =60)
        axes[aa+2].set_title(titles[aa])
        axes[aa+2].set_xlabel("Condition")

        axes[aa].set_title(titles[aa])
        axes[aa].legend(legend)

    plt.savefig('C:/Users/SimonThompson/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/HSV_14Jun_Plots/6_Combined_Whole_Condition_Mean_Perc_Inc_'+ part_size[ab] + '.png')


    # In[334]:


    #7 Combined Means, Well means, percentage increase

    Day_list = ["0619", "0621", "0623", "0625", "0627", "0629", "0701", "0703"]
    df_list = [com_rfp_means, com_gfp_means]
    titles = ["com_rfp", "com_gfp"]
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True,  figsize = (40,10))
    axes = [ax1, ax2, ax3, ax4]
    legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]

    for aa in range(len(df_list)):
        day_jump = 528
        well_jump = 66
        bar = [0]*len(legend)

        for v in range(len(legend)):
            well_j = 66*(v+1)
            c1 = [0]*len(Day_list)
            perc_c1 = [0]*len(Day_list)

            for q in range(len(Day_list)):
                w1 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j-well_jump:well_j-int(((well_jump/3)*2))]["Count"].mean()
                w2 = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)*2):int(well_j-(well_jump/3))]["Count"].mean()
                w3 = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)):well_j]["Count"].mean()

                w_m = [w1, w2, w3]
                c1[q] = np.mean(w_m)
                perc_c1[q] = ((c1[q]-c1[0])/c1[0])*100

                if q == 4:
                    bar[v] = perc_c1[q]

            axes[aa].plot(Day_list, perc_c1)
            axes[aa].set_xlabel("Days")
            axes[aa].set_ylabel("% Increase")

        axes[aa].set_title(titles[aa])
        axes[aa].legend(legend)    

        axes[aa+2].bar(np.arange(len(legend)), bar, color=['blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'grey'])
        axes[aa+2].set_xticks(np.arange(len(legend)))
        axes[aa+2].set_xticklabels(legend, rotation =60)
        axes[aa+2].set_title(titles[aa])
        axes[aa+2].set_xlabel("Condition")

    plt.savefig('C:/Users/SimonThompson/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/HSV_14Jun_Plots/7_Well_Mean_Perc_Inc_'+ part_size[ab] + '.png')


    # In[252]:


    #8 Combined (Li + Triangle), well means, well comparison, percentage increase
    Day_list = ["0619", "0621", "0623", "0625", "0627", "0629", "0701", "0703"]
    df_list = [com_rfp_means, com_gfp_means]
    titles = ["com_rfp"]
    legend = ["US7", "US8", "US9", "US7_8", "US7_9", "US8_9", "US7_8_9", "GFP"]

    fig, ((ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12, ax13, ax14, ax15, ax16)) = plt.subplots(2, 8, sharey=True, sharex=True, figsize = (40,10))
    axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10, ax11, ax12, ax13, ax14, ax15, ax16]

    for aa in range(len(df_list)):
        day_jump = 528
        well_jump = 66

        for v in range(len(legend)):
            well_j = 66*(v+1)

            w1 = [0]*len(Day_list)
            w2 = [0]*len(Day_list)
            w3 = [0]*len(Day_list)
            w_con = [0]*len(Day_list)

            perc_w1 = [0]*len(Day_list)
            perc_w2 = [0]*len(Day_list)
            perc_w3 = [0]*len(Day_list)
            perc_w_con = [0]*len(Day_list)

            w1_std = [0]*len(Day_list)
            w2_std = [0]*len(Day_list)
            w3_std = [0]*len(Day_list)
            w_con_std = [0]*len(Day_list)

            for q in range(len(Day_list)): 
                w1[q] = df_list[aa][q*day_jump:(q+1)*day_jump][well_j-well_jump:well_j-int(((well_jump/3)*2))]["Count"].mean()
                w2[q] = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)*2):int(well_j-(well_jump/3))]["Count"].mean()
                w3[q] = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)):well_j]["Count"].mean()
                w_m = np.mean([w1[q], w2[q], w3[q]])
                w_con[q] = w_m

                perc_w1[q] = ((w1[q]-w1[0])/w1[0])*100
                perc_w2[q] = ((w2[q]-w2[0])/w2[0])*100
                perc_w3[q] = ((w3[q]-w3[0])/w3[0])*100
                perc_w_con[q] = ((w_con[q]-w_con[0])/w_con[0])*100

                w1_std[q] = df_list[aa][q*day_jump:(q+1)*day_jump][well_j-well_jump:well_j-int(((well_jump/3)*2))]["Count"].std()
                w2_std[q] = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)*2):int(well_j-(well_jump/3))]["Count"].std()
                w3_std[q] = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)):well_j]["Count"].std()

            axes[v + (aa*8)].plot(perc_w1, label = 'W1')
            axes[v + (aa*8)].plot(perc_w2, label = 'W2')
            axes[v + (aa*8)].plot(perc_w3, label = 'W3')
            axes[v + (aa*8)].plot(perc_w_con, label = 'Mean')
            axes[v + (aa*8)].legend()
            axes[v + (aa*8)].set_xlabel("Days")

    plt.savefig('C:/Users/SimonThompson/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/HSV_14Jun_Plots/8_Combined_Well_Means_Well_Comp_%Increase_'+ part_size[ab] + '.png')



    # In[234]:


    #9 Li, Well Means, Percentage Increase

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

                w1 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j-well_jump:well_j-int(((well_jump/3)*2))]["Count"].mean()
                w2 = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)*2):int(well_j-(well_jump/3))]["Count"].mean()
                w3 = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)):well_j]["Count"].mean()
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



    plt.savefig('C:/Users/SimonThompson/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/HSV_14Jun_Plots/9_Li_Mean_Perc_Inc_Well'+ part_size[ab] + '.png')


    # In[235]:


    #10, Li and Triangle Means Comparison (Percentage Increase)

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


                w1 = df_list[aa][q*day_jump:(q+1)*day_jump][well_j-well_jump:well_j-int(((well_jump/3)*2))]["Count"].mean()
                w2 = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)*2):int(well_j-(well_jump/3))]["Count"].mean()
                w3 = df_list[aa][q*day_jump:(q+1)*day_jump][int(well_j-(well_jump/3)):well_j]["Count"].mean()

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




    plt.savefig('C:/Users/SimonThompson/Dropbox (UCL - SWC)/Data Files/Simon/Data/HSV Rabies/HSV_HEK_Transfection_14_Jun2018/HSV_14Jun_Plots/10_Li_Triange_Mean_Perc_Inc_Comparison_'+ part_size[ab] + '.png')

        

