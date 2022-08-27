import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt


physics_type = 'dead_oil'  # dead_oil, geothermal
NR_REAL = 10
num_inj_wells = 2
num_prod_wells = 2
use_clean = True
start_plot_time = 1  # First hours are heavily influenced by transient effects and well rates can be a bit misleading
properties_to_plot = ['water rate', 'BHP']  # Note: for geothermal also steam rate, and temp gets plot, for dead_oil also oil rate!

if not use_clean:
    DIR_INPUT = f'ensemble_2_meshes_output_{physics_type}'
else:
    DIR_INPUT = f'ensemble_2_clean_meshes_output_{physics_type}'
BASE_FILENAME = lambda ith_real: f'prod_data_real_{ith_real}.xlsx'


def store_each_well(well_name, num_wells=1, prop='water rate', units='(m3/day)'):
    list_prop = []
    for i in range(1, num_wells + 1):
        list_prop.append(f'{well_name}{i} : {prop} {units}')
    return list_prop

properties = []
properties += store_each_well(well_name='I', num_wells=num_inj_wells, prop='water rate', units='(m3/day)')
properties += store_each_well(well_name='P', num_wells=num_prod_wells, prop='water rate', units='(m3/day)')
properties += store_each_well(well_name='I', num_wells=num_inj_wells, prop='BHP', units='(bar)')
properties += store_each_well(well_name='P', num_wells=num_prod_wells, prop='BHP', units='(bar)')

if physics_type == 'dead_oil':
    properties_to_plot += ['oil rate']
    properties += store_each_well(well_name='I', num_wells=num_inj_wells, prop='oil rate', units='(m3/day)')
    properties += store_each_well(well_name='P', num_wells=num_prod_wells, prop='oil rate', units='(m3/day)')

elif physics_type == 'geothermal':
    properties_to_plot += ['steam rate', 'temperature']
    properties += store_each_well(well_name='I', num_wells=num_inj_wells, prop='steam rate', units='(m3/day)')
    properties += store_each_well(well_name='P', num_wells=num_prod_wells, prop='steam rate', units='(m3/day)')
    properties += store_each_well(well_name='I', num_wells=num_inj_wells, prop='temperature', units='(K)')
    properties += store_each_well(well_name='P', num_wells=num_prod_wells, prop='temperature', units='(K)')

num_props = len(properties)
prod_data_dict = dict()  # BHP Inj, Prod, Rate Inj, Prod,

for i in range(1, NR_REAL + 1):
    loc_data = pd.read_excel(os.path.join(DIR_INPUT, BASE_FILENAME(i)))
    df_loc_data = pd.DataFrame(loc_data)
    loc_prod_data = np.array(loc_data)
    prod_data_dict[i] = np.zeros((loc_prod_data.shape[0], num_props + 1))
    prod_data_dict[i][:, -1] = loc_prod_data[:, df_loc_data.columns.get_loc('time')]
    for j in range(num_props):
        prod_data_dict[i][:, j] = loc_prod_data[:, df_loc_data.columns.get_loc(properties[j])]


font_dict_title = {'family': 'sans-serif',
                   'color': 'black',
                   'weight': 'normal',
                   'size': 14,
                   }

font_dict_axes = {'family': 'monospace',
                  'color': 'black',
                  'weight': 'normal',
                  'size': 14,
                  }

color_list = ['blue', 'red', 'green', 'cyan', 'black', 'yellow', 'purple', 'grey']
for i in range(len(properties_to_plot)):
    loc_property = properties_to_plot[i]
    loc_indices = [id for id, item in enumerate(properties) if loc_property in item]

    fig, axs = plt.subplots(1, 1, figsize=(5, 5), dpi=400, facecolor='w', edgecolor='k')
    for j in range(1, NR_REAL + 1):
        count = 0
        start_id = np.where(prod_data_dict[j][:, -1] > start_plot_time)[0][0]
        for k in range(loc_indices[0], loc_indices[0] + num_inj_wells):
            if j == 1:
                axs.plot(prod_data_dict[j][start_id:, -1], prod_data_dict[j][start_id:, k], color=color_list[count], linewidth=1, label=f'I{count + 1}')
            else:
                axs.plot(prod_data_dict[j][start_id:, -1], prod_data_dict[j][start_id:, k], color=color_list[count], linewidth=1)
            axs.set_xlabel('time [days]', font_dict_axes)
            axs.set_ylabel(loc_property, font_dict_axes)
            axs.set_title(f'Injection wells', fontdict=font_dict_title)
            count += 1
    axs.legend()
    left = 0.05  # the left side of the subplots of the figure
    right = 0.95  # the right side of the subplots of the figure
    bottom = 0.05  # the bottom of the subplots of the figure
    top = 0.95  # the top of the subplots of the figure
    wspace = 0.25  # the amount of width reserved for blank space between subplots
    hspace = 0.25  # the amount of height reserved for white space between subplots
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    axs.tick_params(axis='x', labelsize=16)
    axs.tick_params(axis='y', labelsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(DIR_INPUT, f'Injection_{loc_property}.pdf'))
    plt.show()

    fig, axs = plt.subplots(1, 1, figsize=(5, 5), dpi=400, facecolor='w', edgecolor='k')
    for j in range(1, NR_REAL + 1):
        count = 0
        start_id = np.where(prod_data_dict[j][:, -1] > start_plot_time)[0][0]
        for k in range(loc_indices[0] + num_inj_wells, loc_indices[0] + num_inj_wells + num_prod_wells):
            if j == 1:
                axs.plot(prod_data_dict[j][start_id:, -1], prod_data_dict[j][start_id:, k], color=color_list[count], linewidth=1,
                         label=f'P{count + 1}')
            else:
                axs.plot(prod_data_dict[j][start_id:, -1], prod_data_dict[j][start_id:, k], color=color_list[count], linewidth=1)
            axs.set_xlabel('time [days]', font_dict_axes)
            axs.set_ylabel(loc_property, font_dict_axes)
            axs.set_title(f'Production wells', fontdict=font_dict_title)
            count += 1
    axs.legend()
    left = 0.05  # the left side of the subplots of the figure
    right = 0.95  # the right side of the subplots of the figure
    bottom = 0.05  # the bottom of the subplots of the figure
    top = 0.95  # the top of the subplots of the figure
    wspace = 0.25  # the amount of width reserved for blank space between subplots
    hspace = 0.25  # the amount of height reserved for white space between subplots
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    axs.tick_params(axis='x', labelsize=16)
    axs.tick_params(axis='y', labelsize=16)
    plt.tight_layout()
    plt.savefig(os.path.join(DIR_INPUT, f'Production_{loc_property}.pdf'))
    plt.show()
