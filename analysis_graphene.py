import os
from collections import Counter
import numpy as np
import pandas as pd
import re

path = os.path.dirname(os.path.abspath(__name__))
nomes = [i for i in os.listdir(path) if re.search('results',i)]
nomes


################# DATA IMPORT #################
###############################################
data = {} # all solutions
for nome in nomes:
    data[nome] = pd.read_excel(path+'\\original data\\'+nome).ffill()
    pattern = r'\b\d*' + re.escape('x') + r'\d*\b'
    size = re.findall(pattern, nome, flags=re.IGNORECASE)[0]
    data[nome]['Size'] = size
   
    # # selecting the unique solutions for each size
    # unique_genomes[nome] = data[nome].drop_duplicates(subset='genome', keep='first')
    # unique_genomes[nome] = unique_genomes[nome].sort_values(by='gap')



################### ALL DATA ###################
################################################ 
df = pd.concat(
    [data[nome] for nome in nomes]
)
df.columns = df.columns.str.replace('Unnamed: 0','indrun')
df.columns = df.columns.str.replace('Unnamed: 1','generation')
df['genome'] = list(map(eval,df.genome))
df = df[['Size','indrun', 'generation', 'best_solution', 'genome', 'diag', 'fit', 'gap',
       'ipn', 'homo', 'population']]
################### ALL DATA ###################
################################################




################### CLEANING DATA ###################
#####################################################


sizes = []
indrun = []
generation = []
genome = []
gap = []
ipn = []
homo = []

for k in range(len(df)):
    if len(df.genome.iloc[k]) < 48:
        for j in range(len(df.genome.iloc[k])):
            genome.append(df.genome.iloc[k][j])
            indrun.append(df.indrun.iloc[k])
            generation.append(df.generation.iloc[k])
            sizes.append(df.Size.iloc[k])
            gap.append(eval(df.gap.iloc[k])[j])
            ipn.append(eval(df.ipn.iloc[k])[j])
            homo.append(eval(df.homo.iloc[k])[j])
    else:
        genome.append(df.genome.iloc[k])
        indrun.append(df.indrun.iloc[k])
        generation.append(df.generation.iloc[k])
        sizes.append(df.Size.iloc[k])
        gap.append(df.gap.iloc[k])
        ipn.append(df.ipn.iloc[k])
        homo.append(df.homo.iloc[k])


concentration = list(map(Counter,genome))
hbn = [np.round(i[1]/(i[0] + i[1]),2) for i in concentration]
concentration = [f'$A_{i[0]/(i[0] + i[1]) :.2f} B_{i[1]/(i[0] + i[1]):.2f}$' for i in concentration]
dfs = pd.DataFrame({
    'Size': sizes,
    'Run': indrun,
    'Generation': generation,
    'Concentration': concentration,
    'Band gap [eV]': gap,
    'IPN': ipn,
    'homo': homo,
    'genome': genome,
    'hbn': hbn
})
dfs['chromosome'] = [''.join(list(map(str,i))) for i in dfs.genome]
################### CLEANING DATA ###################
#####################################################




########################################## PLOTS ############################################


################### HISTOGRAM OF SOLUTIONS ###################
##############################################################
import plotly.express as px
import plotly.graph_objects as go

fig = px.histogram(dfs, x="Band gap [eV]", color="Size", marginal="rug", # can be `box`, `violin`
                         hover_data=dfs.columns,
                         text_auto=True,
                        #  labels={'Concentration': [dfs[dfs['Band gap [eV]'] == k]['Concentration'].iloc[0] for k in dfs['Band gap [eV]'].unique()]},
                         color_discrete_map={
                            "6x6": "#1F77B4", # Azul
                            "9x9": "#FF7F0E", # Laranja
                            # "ABC 400": "#2CA02C", # Verde
                            # "ABC 500": "#D62728", # Vermelho
                            }
                        )

fig.update_layout(
# height=800,
# width=1300,
template='plotly_white',
font_family="Arial",
font_color="black",
title_font_family="Arial",
# title_font_color="red",
legend_title= "Size",
font_size = 18,
# title = 'History of evolution',
# xaxis_title="Run",
yaxis_title="Number of solutions",
xaxis = dict(
    # tickmode = 'linear',
    mirror=True,
    ticks='outside',
    showline=True,
    tickformat=".1f"
    ),
yaxis = dict(
    # tickmode = 'linear',
    # tick0 = 22,
    # dtick = 10,
    mirror=True,
    ticks='outside',
    showline=True,
    
    ),

legend=dict(
    # orientation="h",
    entrywidth=130,
    # yanchor="middle",
    # y=1.02,
    # xanchor="left",
    # x=0.9,
    font_size = 16
    )
,

)


fig.update_traces(opacity=0.85)
fig.update_layout(
    barmode='stack',
    # barmode='overlay',
    # bargap=0.2, # gap between bars of adjacent location coordinates
    # bargroupgap=0.1 # gap between bars of the same location coordinates
    )

fig.show()

################### HISTOGRAM OF SOLUTIONS ###################
##############################################################
histogram = fig


################### HEATMAP DENSITY ###################
#######################################################

import plotly.express as px
import plotly.graph_objects as go


density = dfs[dfs['Band gap [eV]'] <= 1.0][dfs['Band gap [eV]'] >= 0.5]
density.sort_values(by=['Band gap [eV]', 'hbn'])

fig = px.density_heatmap(
    density,
    x="Band gap [eV]",
    y="hbn",
    # z = 'Size',
    # marginal_x="histogram",
    # marginal_y="histogram",
    # range_x=[0.5,1.0],
    # histfunc= 'count',
    labels = dict(z ="Productivity"),
    text_auto=True
    )


fig.update_layout(
# height=800,
# width=1300,
template='plotly_white',
font_family="Arial",
font_color="black",
title_font_family="Arial",
# title_font_color="red",
# legend_title= "Size",
font_size = 16,
# title = 'History of evolution',
# xaxis_title="Run",
yaxis_title="h-BN Concentration",
xaxis = dict(
    mirror=True,
    ticks='outside',
    showline=True,
    # tickformat=".2s"
    ),
yaxis = dict(
    # tickmode = 'linear',
    mirror=True,
    ticks='outside',
    showline=True,
    tickformat = 'p'
    ),
coloraxis_colorbar_title_text = 'Number of solutions',
legend=dict(
    # orientation="h",
    entrywidth=130,
    # yanchor="middle",
    # y=1.02,
    # xanchor="left",
    # x=0.9,
    font_size = 14
    )
,
)


fig.show()

################### HEATMAP DENSITY ###################
#######################################################
heatmap = fig



dfunique = dfs.drop_duplicates(subset=['chromosome'])



###################### HEATMAP DENSITY ####################
######################### UNIQUES #########################

import plotly.express as px
import plotly.graph_objects as go


density = dfunique[dfunique['Band gap [eV]'] <= 1.0][dfunique['Band gap [eV]'] >= 0.5]
density.sort_values(by=['Band gap [eV]', 'hbn'])

fig = px.density_heatmap(
    density,
    x="Band gap [eV]",
    y="hbn",
    # z = 'Size',
    # marginal_x="histogram",
    # marginal_y="histogram",
    # range_x=[0.5,1.0],
    # histfunc= 'count',
    labels = dict(z ="Productivity"),
    text_auto=True
    )


fig.update_layout(
# height=800,
# width=1300,
template='plotly_white',
font_family="Arial",
font_color="black",
title_font_family="Arial",
# title_font_color="red",
# legend_title= "Size",
font_size = 16,
# title = 'History of evolution',
# xaxis_title="Run",
yaxis_title="h-BN Concentration",
xaxis = dict(
    mirror=True,
    ticks='outside',
    showline=True,
    # tickformat=".2s"
    ),
yaxis = dict(
    # tickmode = 'linear',
    mirror=True,
    ticks='outside',
    showline=True,
    tickformat = 'p'
    ),
coloraxis_colorbar_title_text = 'Number of solutions',
legend=dict(
    # orientation="h",
    entrywidth=130,
    # yanchor="middle",
    # y=1.02,
    # xanchor="left",
    # x=0.9,
    font_size = 14
    )
,
)

# fig.add_trace(go.Scatter(
#     x=density[density['Size'] == '6x6']['Band gap [eV]'],
#     y=density[density['Size'] == '6x6']['hbn'],
#     mode='markers',
#     showlegend=False,
#     marker=dict(
#         symbol='x',
#         opacity=0.7,
#         color='white',
#         size=8,
#         line=dict(width=1),
#     )
# ))
# fig.add_trace(go.Scatter(
#     x=density[density['Size'] == '9x9']['Band gap [eV]'],
#     y=density[density['Size'] == '9x9']['hbn'],
#     mode='markers',
#     showlegend=False,
#     marker=dict(
#         symbol='circle',
#         opacity=0.7,
#         color='white',
#         size=8,
#         line=dict(width=1),
#     )
# ))

fig.show()

###################### HEATMAP DENSITY ####################
######################### UNIQUES #########################
heatmapunique = fig


gf = density[density['Band gap [eV]'] <= 0.975][density['Band gap [eV]'] >= 0.875][density['hbn'] >= 0.48][density['hbn'] <= 0.52].sort_values(by=['IPN'])

x6 = gf[gf.Size == '6x6'].genome.iloc[0]
x9 = gf[gf.Size == '9x9'].genome.iloc[0]


##### calculations #####

import peakutils as pk
from graphene_unitcell import hmo, tb_matrix_unit

energies6 = hmo(tb_matrix_unit(x6))
energies9 = hmo(tb_matrix_unit(x9))



def density_of_states(energy_levels):
    
    # Calculate the DOS by counting the energy levels
    dos = {}
    for energy in energy_levels:
        if energy in dos:
            dos[energy] += 1
        else:
            dos[energy] = 1
            
    # Extract unique energy values and corresponding counts
    unique_energy_values = list(dos.keys())
    states = [dos[energy] for energy in unique_energy_values]
    
    return unique_energy_values, states


def gaussian_peaks(centers, states, sigma,
                   x0 = None, xf = None, delta = None):
    
    if delta == None:
        x = np.linspace(x0, xf, num=100000) # type: ignore
    else:
        x = np.linspace(np.min(centers)-delta,np.max(centers)+delta, 100000)
        
    # centers = unique_energy_values
    y = 0

    for i in range(len(centers)):
        y += pk.gaussian(x, states[i], centers[i], sigma)
        
    return x,y

x1, y1 = density_of_states(energies6['energies'])
x1, y1 = gaussian_peaks(centers = x1, states = y1, sigma = 0.02, x0=-3, xf=1.5)

x2, y2 = density_of_states(energies9['energies'])
x2, y2 = gaussian_peaks(centers = x2, states = y2, sigma= 0.02,x0=-3,xf=1.5)





##### figures #####

import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

matplotlib.rc('xtick', labelsize=10) 
matplotlib.rc('ytick', labelsize=10) 
plt.figure(figsize=(3.54,2.188976), dpi = 600)

fig, axes = plt.subplots(nrows=2, sharex=True)
fig.subplots_adjust(hspace=0)

############ SUBPLOT 1 ############
###################################
ax1 = axes[0]
# divider1 = make_axes_locatable(ax1) # divide a subplot
# ax12 = divider1.new_vertical(size='100%',pad=0.12)
# fig.add_axes(ax12)

# hide the spines between ax@ and ax@2
ax1.plot(x1,y1)
# ax1.spines.top.set_visible(False)
# ax12.plot(x1,y1)
# ax12.spines.bottom.set_visible(False)
# ax12.tick_params(labelbottom=False, bottom = False)

# zoom-in / limit the view to different portions of the data
ax1.set_ylim(0,3)
# ax12.set_ylim(100,230)


# add legends and lines
ax1.axvline(energies6['energies'][energies6['homo']],linewidth=2, color='r',linestyle='--', label = 'HOMO')
ax1.axvline(energies6['energies'][energies6['homo']+1],linewidth=2, color='b',linestyle='--', label = 'LUMO')
# ax12.axvline(energies6['energies'][energies6['homo']],linewidth=2, color='r',linestyle='--', label = 'HOMO')
# ax12.axvline(energies6['energies'][energies6['homo']+1],linewidth=2, color='b',linestyle='--', label = 'LUMO')
ax1.legend(fontsize = 12, loc = 'upper left')
# ############ SUBPLOT 1 ############
###################################


############ SUBPLOT 2 ############
###################################
ax2 = axes[1]
# divider2 = make_axes_locatable(ax2) # divide a subplot
# ax22 = divider2.new_vertical(size='100%',pad=0.12)
# fig.add_axes(ax22)

# hide the spines between ax@ and ax@2
ax2.plot(x2,y2)
# ax2.spines.top.set_visible(False)
# ax22.plot(x2,y2)
# ax22.spines.bottom.set_visible(False)
# ax22.tick_params(labelbottom=False, bottom = False)

# zoom-in / limit the view to different portions of the data
ax2.set_ylim(0,4)
# ax22.set_ylim(100,200)


# add legends and lines
ax2.axvline(energies9['energies'][energies9['homo']],linewidth=2, color='r',linestyle='--', label = 'HOMO')
ax2.axvline(energies9['energies'][energies9['homo']+1],linewidth=2, color='b',linestyle='--', label = 'LUMO')

# ax22.axvline(energies9['energies'][energies9['homo']],linewidth=2, color='r',linestyle='--', label = 'HOMO')
# ax22.axvline(energies9['energies'][energies9['homo']+1],linewidth=2, color='b',linestyle='--', label = 'LUMO')
# ax2.legend(fontsize = 12, loc = 'bottom left')
############ SUBPLOT 2 ############
###################################


ax1.grid()
ax2.grid()

ax1.set_yticks([])
ax2.set_yticks([])


# Set common labels
fig.text(0.5, 0.04, 'Energy [\u03B2]', ha='center', va='center', fontsize= 15)
fig.text(0.05, 0.5, 'Density of States', ha='center', va='center', rotation='vertical',fontsize= 15)
dos = fig





sf = dfunique[dfunique['Band gap [eV]'] <= 1.0][dfunique['Band gap [eV]'] >= 0.5].sort_values(by=['Size', 'IPN'])



if __name__ == '__main__':
    dos.savefig('density of states graphene.svg', format = 'svg', bbox_inches='tight',dpi=600)
    sf.to_excel('unique solution within range, sorted by IPN.xlsx')
    histogram.write_image("histogram graphene.svg",width=720, height=400, scale = 1)
    heatmap.write_image("heatmap graphene.svg",width=720, height=400, scale = 1)