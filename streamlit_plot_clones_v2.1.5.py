import itertools

import streamlit
import streamlit.components.v1 as streamlit_components_v1
import pandas
import os
import plotly_express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.tools as plotly_tools
import json
import mpld3

import matplotlib.pyplot as matplotlib_pyplot
from matplotlib.colors import ListedColormap
from pyfish import fish_plot, process_data, setup_figure

script_dir = os.path.dirname(__file__)

__version__ = '2.1.6'

streamlit.set_page_config(layout="wide")

# User input - File upload - starting point
uploaded_file = streamlit.sidebar.file_uploader("Choose a file to start", key='file_uploader')


def hex_2_rgb(hex_string):
    """
    Returns an RGB color as a list, so an alpha value can be appended
    :param hex_string:
    :return:
    """
    r_hex = hex_string[1:3]
    g_hex = hex_string[3:5]
    b_hex = hex_string[5:7]
    return [int(r_hex, 16), int(g_hex, 16), int(b_hex, 16)]  # We return a list, so we can append


def export_to_svg_ui(plotly_figure):
    svg_export_ui = streamlit.sidebar.expander('SVG Export')
    name = svg_export_ui.text_input('Custom file name')
    file_name = 'svg_file.svg'
    if name:
        if name.endswith('.svg'):
            file_name = name
        else:
            file_name = name + '.svg'

    svg_export_ui.markdown(f'File will be saved as:\n:red[{file_name}]')

    svg_export_ui.download_button(
        label="Download SVG",
        data=plotly_figure.to_image(format='svg', validate=True),
        file_name=file_name,
        mime="text/svg",
        icon=":material/download:",
    )


def get_plot_theme_and_frame_size(container):
    plot_theme = container.selectbox('Plot theme', ['plotly', 'plotly_white', 'plotly_dark', 'ggplot2', 'seaborn', 'simple_white', 'presentation', 'xgridoff', 'ygridoff', 'gridon', 'none'])

    frame_layout = streamlit.sidebar.expander('Frame layout')
    frame_layout.write('If any of the plots are "cut off", you can increase the size of the plot display frame here - rarely necessary.')
    frame_width = frame_layout.slider('Frame width', min_value=1000, max_value=5000, step=100, value=3000)
    frame_height = frame_layout.slider('Frame height', min_value=500, max_value=2500, step=100, value=900)

    return plot_theme, frame_width, frame_height


if uploaded_file is not None:
    # streamlit.write(streamlit.session_state)

    patient_name = os.path.splitext(os.path.basename(uploaded_file.name))[0]
    data_frame = pandas.read_csv(uploaded_file, sep='\t')

    # Sort data frame relative to cluster names (in order to get legends sorted as they will be in the order they are added)
    # data_frame = data_frame.sort_values(by='Cluster', key=lambda col: col.map(lambda x: int(x.split('_')[1])))

    column_names = data_frame.columns
    vaf_columns = [name for name in column_names if (name.startswith('VAF') and not name.startswith('VAF_CCF'))]
    sample_names = [name.split('_')[1] for name in vaf_columns]

    all_clusters = data_frame.Cluster.unique()
    clusters = all_clusters

    # Combine sample names

    # Since we do not know the interrelation of samples - i.e. which is primary tumor and which is metastasis etc. we make
    # no attempt to organize combinations in terms of which sample name is first/last in combinations. As this, however,
    # will effect what axis each sample goes to, we implement the ability to switch axes on the plots later in the script.

    #################################################################################
    #                              Select plot type                                 #
    #################################################################################
    # Figure out if Fish plot should be an option - if fish plot data exist in the dataframe
    plot_types = ['Dot plot', 'Heat map', '3D line plot', '3D surface plot', 'ctDNA plot']
    try:
        tree_branches = [x for x in data_frame.tree_branches.unique() if isinstance(x, str)][0].split(';')
        clone_proportions = json.loads([x for x in data_frame.clone_proportions.unique() if isinstance(x, str)][0])
        plot_types.append('Fish plot')
    except AttributeError as error:
        if "tree_branches" not in str(error):
            raise error

    plot_type_ui = streamlit.sidebar.expander('Plot type')
    plot_type = plot_type_ui.radio('Select plot type', plot_types, key='plot_type')

    if not plot_type == 'Fish plot':
        #################################################################################
        #                           SELECT CLUSTERS EXPANDER                            #
        #################################################################################

        # Define UI

        def check_all():
            for cluster in all_clusters:
                streamlit.session_state[cluster] = True


        def uncheck_all():
            for cluster in all_clusters:
                streamlit.session_state[cluster] = False


        select_clusters = streamlit.sidebar.expander('Select clusters')
        check_all_button = select_clusters.button('Check all', 'check_all_button', on_click=check_all)
        uncheck_all_button = select_clusters.button('Uncheck all', 'uncheck_all_button', on_click=uncheck_all)

        show_cluster = {}
        for cluster in all_clusters:
            if cluster not in streamlit.session_state:
                streamlit.session_state[cluster] = True
            show_cluster[cluster] = select_clusters.checkbox(cluster, key=cluster)

        # Filter data frame based on show_cluster values
        clusters_to_include = [cluster for cluster in show_cluster if show_cluster[cluster]]
        data_frame = data_frame[data_frame['Cluster'].isin(clusters_to_include)]

        clusters = clusters_to_include
    else:
        clusters = all_clusters

    cluster_colors = list(data_frame.Color.unique())

    if plot_type == 'Dot plot':

        pairwise_sample_combinations = itertools.combinations(sample_names, 2)

        display_combinations = {'{}_vs_{}'.format(combination[0], combination[1]): list(combination) for combination in
                                pairwise_sample_combinations}
        number_of_plots = len(display_combinations)

        if number_of_plots > 1:

            # Add the choice of multiplot to the samples dropdown
            display_combinations['MultiPlot'] = 'MultiPlot'

            # User input - Which samples to plot
            sample_combination = display_combinations[
                streamlit.sidebar.selectbox('Please select which samples to compare', display_combinations.keys())]
        else:
            sample_combination = list(display_combinations.values())[0]

        # Data filtering
        # data_filtering = streamlit.sidebar.expander('Data filters')
        #
        # # User input
        # min_vaf = data_filtering.slider('Minimal MAF', min_value=0.0, max_value=1.0, value=0.0, step=0.01)

        # User input - Which data type to plot
        data_type = streamlit.sidebar.radio('Select data type', ('VAF', 'cluster_CCF', 'VAF_CCF'))

        #################################################################################
        #                           VISUAL APPEARANCE EXPANDER                          #
        #################################################################################

        # Define Expander for visual appearance
        visual_appearance = streamlit.sidebar.expander('Plot - visual appearance')

        plot_theme, frame_width, frame_height = get_plot_theme_and_frame_size(visual_appearance)

        # User input - dot size
        dot_size = visual_appearance.selectbox('Dot size...', range(5, 21), index=3)

        # User input - display/hide dot periphery line
        display_dot_periphery_line = visual_appearance.checkbox('Toggle dot edge-lines', value=False)

        if display_dot_periphery_line:
            marker = {
                'size': dot_size,
                'line': {
                    'width': 2,
                    'color': 'DarkSlateGrey'
                }
            }
        else:
            marker = {
                'size': dot_size,
            }

        # Define ranges for axes (toggle between range for VAF plot or CCF plot)
        x_y_ranges = ([-0.05, 1.05], [-0.1, 2.1])
        range_x = x_y_ranges[0] if data_type == 'VAF' else x_y_ranges[1]
        range_y = x_y_ranges[0] if data_type == 'VAF' else x_y_ranges[1]

        # Multiplot
        if sample_combination == 'MultiPlot':

            # User input
            grid_columns = visual_appearance.slider('Number of grid columns', min_value=2, max_value=8, value=3, step=1)

            # User input
            inter_space = visual_appearance.slider('Space between plots', min_value=0.05, max_value=0.5, value=0.2,
                                                   step=0.05)

            grid_rows = int(number_of_plots / grid_columns) + 1 if number_of_plots % grid_columns else int(
                number_of_plots / grid_columns)

            subplot_width = visual_appearance.number_input('Sub-plot width', min_value=100, max_value=1000, value=400, step=100)

            plot_height_subtraction = visual_appearance.number_input('Plot height subtraction', min_value=-200, max_value=500, value=50, step=1)

            subplot_titles = list(display_combinations.keys())[:-1]

            # Define plot settings for axis flipping
            plot_settings = streamlit.sidebar.expander('Flip axes')
            flip_combination_axes = {}
            for title in subplot_titles:
                flip_combination_axes[title] = plot_settings.checkbox(title, value=False, key=title)

            for sample_combination in flip_combination_axes:
                if flip_combination_axes[sample_combination]:
                    display_combinations[sample_combination] = display_combinations[sample_combination][::-1]
            h_space = inter_space / grid_columns
            v_space = inter_space / grid_rows

            figure = make_subplots(cols=grid_columns, rows=grid_rows, subplot_titles=subplot_titles,
                                   horizontal_spacing=h_space, vertical_spacing=v_space)

            row = 1
            col = 1
            plot_number = 1
            legend_groups = set()
            for sample_combination in list(display_combinations.values())[:-1]:

                # We currently have no structure to guarantee primary tumors on x axes and metastases on y since we
                x_y_axes = (data_type + '_' + sample_combination[0], data_type + '_' + sample_combination[1])

                px_figure = px.scatter(
                    data_frame,
                    x=x_y_axes[0],
                    y=x_y_axes[1],
                    range_x=range_x,
                    range_y=range_y,
                    color="Cluster",
                    color_discrete_sequence=cluster_colors,
                    facet_col="Cluster",
                    hover_data={
                        # x_y_axes[0]: False,  # Displays VAF value
                        # x_y_axes[1]: False,
                        'ref_counts_' + sample_combination[0]: True,
                        'alt_counts_' + sample_combination[0]: True,
                        'ref_counts_' + sample_combination[1]: True,
                        'alt_counts_' + sample_combination[1]: True,
                        'major_cn_' + sample_combination[0]: True,
                        'minor_cn_' + sample_combination[0]: True,
                        'major_cn_' + sample_combination[1]: True,
                        'minor_cn_' + sample_combination[1]: True,
                        'Cluster': True,
                        "Mutation": True,
                        "Variant_Type": True,
                        'Impact': True,
                        'Gene': True,
                    },
                )

                for trace in px_figure['data']:

                    # Avoid having the legend duplicated for each added subplot
                    if trace['legendgroup'] in legend_groups:
                        trace['showlegend'] = False
                    else:
                        legend_groups.add(trace['legendgroup'])
                    figure.add_trace(trace, row=row, col=col)

                plot_number += 1
                row = int(plot_number / grid_columns) + 1 if plot_number % grid_columns else int(plot_number / grid_columns)
                col = plot_number % grid_columns or grid_columns

            figure.update_traces(
                marker=marker,
                # selector=dict(mode='markers'),

            ).update_layout(
                hoverlabel_align='left',  # Necessary for streamlit to make text for all labels align left
                width=subplot_width * grid_columns + (grid_columns - 1) * h_space,
                height=(subplot_width - plot_height_subtraction) * grid_rows + (grid_rows - 1) * v_space,
                template=plot_theme
            ).update_yaxes(
                range=range_y
            ).update_xaxes(
                range=range_x
            )

        # Single Plot
        else:
            # User input - plot width
            plot_width = visual_appearance.number_input('Plot width', min_value=200, max_value=2000, value=700, step=50)

            plot_height_subtraction = visual_appearance.number_input('Plot height subtraction', min_value=0, max_value=1000, value=108, step=1)

            # User input - Flip axes
            flip_axes = streamlit.checkbox('Flip axes', value=False)
            orig_sample_combination = sample_combination
            if flip_axes:
                sample_combination = sample_combination[::-1]
            else:
                sample_combination = orig_sample_combination

            # If for some reason variants are uninformative in some samples (var + ref_counts = 0), we may indicate that by
            # finding indexes for those cases and individually styling each corresponding marker in the trace, using:
            # https://stackoverflow.com/questions/70275607/how-to-highlight-a-single-data-point-on-a-scatter-plot-using-plotly-express
            #
            # # User input - indicate uninformative VAFs
            # indicate_uninformative_vafs = streamlit.checkbox('Indicate potentially uninformative VAFs', value=False)
            #
            # if indicate_uninformative_vafs:
            #     uninformative_vaf_indexes = data_frame.index[data_frame['ref_counts_' + sample_combination[0]] == 0].tolist()
            #     # uninformative_vaf_indexes = data_frame.index[
            #     #     (
            #     #         (data_frame['ref_counts_' + sample_combination[0]] == 0) &
            #     #         (data_frame['alt_counts_' + sample_combination[0]]) == 0
            #     #     ) |
            #     #     (
            #     #         (data_frame['ref_counts_' + sample_combination[1]] == 0) &
            #     #         (data_frame['alt_counts_' + sample_combination[1]] == 0)
            #     #     )
            #     # ].tolist()
            #
            #     streamlit.write('ref_counts_' + sample_combination[0])
            #     streamlit.write(uninformative_vaf_indexes)
            #
            # However, this does not seem to be a problem with current datasets, so not implemented

            # We currently have no structure to guarantee primary tumors on x axes and metastases on y. We could implement a manual axis flip
            x_y_axes = (data_type + '_' + sample_combination[0], data_type + '_' + sample_combination[1])

            figure = px.scatter(data_frame, x=x_y_axes[0], y=x_y_axes[1],
                                range_x=range_x,
                                range_y=range_y,
                                color="Cluster",
                                color_discrete_sequence=cluster_colors,
                                width=plot_width,
                                height=plot_width - plot_height_subtraction,
                                # facet_col='Sample',
                                hover_data={
                                    # x_y_axes[0]: False,
                                    # x_y_axes[1]: False,
                                    'ref_counts_' + sample_combination[0]: True,
                                    'alt_counts_' + sample_combination[0]: True,
                                    'ref_counts_' + sample_combination[1]: True,
                                    'alt_counts_' + sample_combination[1]: True,
                                    'major_cn_' + sample_combination[0]: True,
                                    'minor_cn_' + sample_combination[0]: True,
                                    'major_cn_' + sample_combination[1]: True,
                                    'minor_cn_' + sample_combination[1]: True,
                                    'Cluster': True,
                                    "Mutation": True,
                                    "Variant_Type": True,
                                    'Impact': True,
                                    'Gene': True,
                                },
                                )

            figure.update_traces(
                marker=marker,
                selector=dict(mode='markers'),
            ).update_layout(
                hoverlabel_align='left',  # Necessary for streamlit to make text for all labels align left
                template=plot_theme
            )

        # svg_data = figure.to_image(format='svg')

        # streamlit.plotly_chart(figure, theme=plot_theme, use_container_width=False) # Old way
        streamlit_components_v1.html(figure.to_html(), width=frame_width, height=frame_height, scrolling=True)

        export_to_svg_ui(figure)

    elif plot_type == 'Heat map':
        # Add UI for manipulating the heatmap
        visual_appearance = streamlit.sidebar.expander('Plot - visual appearance')

        plot_theme, frame_width, frame_height = get_plot_theme_and_frame_size(visual_appearance)

        #################################################################################
        #                               Set gap between tiles                           #
        #################################################################################

        tile_gap = visual_appearance.number_input('Gap between tiles', min_value=0, max_value=4, value=1, step=1)


        #################################################################################
        #                              Set plot height                                 #
        #################################################################################
        number_of_rows = len(data_frame)

        def change_plot_height():
            streamlit.session_state.user_set_plot_height = streamlit.session_state.user_set_tile_height * number_of_rows

        def change_row_height():
            streamlit.session_state.user_set_tile_height = round(streamlit.session_state.user_set_plot_height / number_of_rows)


        # Set default values useing the session-state API rather than the value argument
        if 'user_set_tile_height' not in streamlit.session_state:
            streamlit.session_state.user_set_tile_height = 18
        if 'user_set_plot_height' not in streamlit.session_state:
            streamlit.session_state.user_set_plot_height = 18 * number_of_rows

        user_set_tile_height = visual_appearance.number_input('Tile height', min_value=1, max_value=25, step=1, key='user_set_tile_height', on_change=change_plot_height)

        user_set_plot_height = visual_appearance.number_input('Plot height', key='user_set_plot_height', on_change=change_row_height)

        user_set_plot_width = visual_appearance.number_input('Plot width', min_value=200, max_value=2000, value=800, step=10)
        # plot_height = user_set_plot_height
        #
        # visual_appearance.write('Current plot height: ')
        # visual_appearance.write(plot_height)

        # Strategy: use Heatmap from plotly.graph_objects - one trace per clone and potentially zero space between plots

        #################################################################################
        #                              Space between color bars                         #
        #################################################################################

        color_bars_top_pos = visual_appearance.number_input('Color bars top position', min_value=0.9, max_value=1.1, value=0.95, step=0.001, format="%0.3f")

        color_bars_distance_factor = visual_appearance.number_input('Color bars vertical spacing factor', min_value=0.001, max_value=0.5, value=0.01, step=0.001, format="%0.3f")

        clusters = data_frame.Cluster.unique()

        genes_per_cluster = {}

        for cluster in clusters:
            genes_per_cluster[cluster] = len(data_frame[data_frame['Cluster'] == cluster])

        subplot_heights = [genes_per_cluster[cluster] / number_of_rows for cluster in clusters]  # In fractions of the whole plot

        # For vertical bars (bad if one or more clusters have few variants):
        # Placing color bars at approximately the middle of each subplot (via fractions of the whole plot)
        # New traces are added top down, so we need to reflect that on the factor - i.e. we start high and go to low in the y-position
        # color_bars_y_position = [sum(subplot_heights[i + 1:]) + subplot_heights[i] / 2 for i in range(len(subplot_heights))]

        color_bars_y_position = [color_bars_top_pos - i * color_bars_distance_factor for i in range(len(clusters))]

        color_bars_x_position = visual_appearance.number_input('Color bars horizontal position', min_value=1.01, max_value=1.5, value=1.12, step=0.01)

        figure_right_margin = visual_appearance.slider('Figure right margin', value=120, min_value=0, max_value=1000, step=20)


        #################################################################################
        #                               Remove gene names                               #
        #################################################################################

        remove_gene_names = visual_appearance.checkbox('Remove gene names', value=False)

        figure = make_subplots(
            len(clusters),
            1,
            shared_xaxes=True,
            row_heights=subplot_heights,
            vertical_spacing=0,
            x_title='Sample',
            y_title='Gene',
        )
        figure.layout.annotations[0].yshift = -30
        figure.layout.annotations[0].text = '<b>Sample</b>'
        figure.layout.annotations[1].xshift = -60
        figure.layout.annotations[1].text = '<b>Gene</b>'

        figure.update_yaxes(
            showgrid=False,
            # title_font_size=20,
        )
        figure.update_xaxes(
            showgrid=False,
            tickcolor='rgba(0,0,0,0)',
        )

        # The idea is to display rgba values using the cluster colors in turn as colors and the VAF as a transparency.

        color_values = []

        vaf_columns = []
        for sample in sample_names:
            vaf_columns.append('VAF_' + sample)

        # Because each subplot gets its own yaxis config (see: https://stackoverflow.com/questions/72464495/plotly-how-to-change-ticks-in-a-subplot), we update nticks for each (we set it to the row number of the full figure, it seems the actual number is set to the number for the particular cluster)
        yaxis_nticks_layout_dict = {}

        for i, cluster in enumerate(clusters):
            vaf_values = []

            # Store the rgb color for this cluster
            cluster_rgb = hex_2_rgb(cluster_colors[i])

            # Get dataframe for just this cluster
            cluster_df = data_frame[data_frame['Cluster'] == cluster]
            cluster_y_axis_genes = cluster_df.Gene.tolist()

            yaxis_nticks_layout_dict['yaxis{}_nticks'.format('' if i == 0 else str(i + 1))] = len(cluster_y_axis_genes)

            if remove_gene_names:
                yaxis_nticks_layout_dict['yaxis{}_showticklabels'.format('' if i == 0 else str(i + 1))] = False
                yaxis_nticks_layout_dict['yaxis{}_tickcolor'.format('' if i == 0 else str(i + 1))] = 'rgba(0,0,0,0)'

            for index, row in cluster_df.iterrows():
                row_vaf_values = []
                for vaf_column in vaf_columns:
                    row_vaf_values.append(row[vaf_column])
                vaf_values.append(row_vaf_values)

            color_bar = dict(
                title=dict(
                    text='_'.join(cluster.split('_')[:2]),
                    side='top'
                ),
                tickmode='auto',
                lenmode='pixels',
                thicknessmode='pixels',
                len=70,
                thickness=10,
                y=color_bars_y_position[i],
                x=color_bars_x_position,
                orientation='h'  #  Makes each color bar horizontal on the otherwise vertical legend
            )

            figure.add_trace(
                go.Heatmap(
                    x=sample_names,
                    y=cluster_y_axis_genes,
                    z=vaf_values,
                    colorscale=[
                        [0, 'rgba({},{},{},0)'.format(cluster_rgb[0], cluster_rgb[1], cluster_rgb[2])],
                        [1, 'rgba({},{},{},1)'.format(cluster_rgb[0], cluster_rgb[1], cluster_rgb[2])]
                    ],
                    colorbar=color_bar,
                    # showscale=False,
                    xgap=tile_gap,  # Gap between heatmap tiles
                    ygap=tile_gap,
                ),
                row=(i + 1), col=1,
            )


        # height_factor = streamlit.slider('Plot height', value=18, min_value=1, max_value=25, step=1)

        figure.update_layout(
            # plot_bgcolor='white',
            height=user_set_plot_height or number_of_rows * user_set_tile_height,
            width=user_set_plot_width,
            template=plot_theme,
            margin=dict(r=figure_right_margin, t=0, l=0, b=0)
        )

        # Also updating the layout for the y axes ticks
        figure.update_layout(yaxis_nticks_layout_dict)

        # figure.layout.yaxis.showticklabels = False

        # streamlit.plotly_chart(figure, use_container_width=False, sharing="streamlit", theme=None)
        streamlit_components_v1.html(figure.to_html(), width=frame_width, height=frame_height, scrolling=True)

        export_to_svg_ui(figure)

    elif plot_type == '3D line plot':
        # First we need to reformat the data. We create a new data frame with four columns: Sample, Gene, VAF, Cluster
        new_data = {
            'Sample': [],
            'Gene': [],
            'VAF': [],
            'Cluster': [],
        }

        gene_column = data_frame.Gene.tolist()
        for vaf_column_name in vaf_columns:
            sample_name = vaf_column_name.split('_')[1]
            vaf_column = data_frame[vaf_column_name].tolist()
            gene_column = list(range(len(gene_column)))
            cluster_column = data_frame.Cluster.tolist()
            sample_column = len(vaf_column) * [sample_name]
            new_data['Sample'] += sample_column
            new_data['Gene'] += gene_column
            new_data['VAF'] += vaf_column
            new_data['Cluster'] += cluster_column

        # We replace the gene names with numbers - to account for genes that are represented more than once
        # tick_text = gene_column

        new_df = pandas.DataFrame(new_data)

        # We sort relative to cluster
        sorted_new_df = new_df.sort_values(
            by='Cluster',
            kind='mergesort'
        )

        sample_frames = {}
        for sample in sample_names:
            sample_frames[sample] = sorted_new_df[sorted_new_df['Sample'].isin([sample])]

        # Add the first figure in order to include the trace
        # plotly_figures = [
        #     px.line_3d(sample_frames[sample], x='Sample', y='Gene', z='VAF', color='Cluster', template='plotly')
        # ]
        plotly_figures = []
        for sample in list(sample_frames.keys()):
            plotly_figures.append(px.line_3d(sample_frames[sample], x='Sample', y='Gene', z='VAF', color='Cluster', color_discrete_sequence=cluster_colors))

        go_data = plotly_figures[0].data

        for plotly_figure in plotly_figures[1:]:
            go_data += plotly_figure.data

        figure = go.Figure(data=go_data)

        # Avoid having cluster legends repeated for each trace
        legend_groups = set()
        for trace in figure['data']:

            # Avoid having the legend duplicated for each added subplot
            if trace['legendgroup'] in legend_groups:
                trace['showlegend'] = False
            else:
                legend_groups.add(trace['legendgroup'])

        figure.update_layout(
            width=1000,
            height=1000,
            scene=dict(
                xaxis=dict(
                    title_text='Sample'
                ),
                yaxis=dict(
                    ticktext=gene_column,
                    title_text='"Variant"'
                    # nticks=len(gene_column)
                ),
                zaxis=dict(
                    title_text='VAF'
                ),
            ),
            template='plotly'
        )
        streamlit.plotly_chart(figure, theme='streamlit', use_container_width=False)
    elif plot_type == '3D surface plot':
        from scipy.signal import savgol_filter

        #################################################################################
        #                                 Edit plot UI                                  #
        #################################################################################
        edit_plot_ui = streamlit.sidebar.expander('Edit plot')

        #################################################################################
        #                             Manipulate plot UI                                #
        #################################################################################

        manipulate_plot = edit_plot_ui.radio('Manipulate plot', ['Original', 'Savitzky-Golay smoothing', 'Peak sorting'])

        # Separate dataframe into dataframes for each cluster. Since we want each cluster to take up a separate space at
        # the X/Y dimensions, we must for each cluster fill the other "clusters" spaces with NaN values.
        if manipulate_plot == 'Savitzky-Golay smoothing':

            #################################################################################
            #                  Set window size for Savitzky-Golay smoothing                 #
            #################################################################################

            window_size = edit_plot_ui.slider('Window size', min_value=2, max_value=150, value=5)

        elif manipulate_plot == 'Peak sorting':
            sort_type = edit_plot_ui.radio('Sort by sample having the', ['highest VAF peak', 'highest overall VAF'])

        cluster_dataframes = {}
        for i, cluster in enumerate(clusters):

            # get VAF columns data for current cluster only
            cluster_data = data_frame[data_frame['Cluster'].isin([cluster])].loc[:, vaf_columns]
            if manipulate_plot == 'Savitzky-Golay smoothing':
                if len(cluster_data) > 1:
                    cluster_data = cluster_data.apply(lambda x: savgol_filter(x, min(window_size, len(cluster_data)),1), axis=0)
            elif manipulate_plot == 'Peak sorting':
                if sort_type == 'highest VAF peak':
                    max_vaf_sample = cluster_data.max().sort_values(ascending=False).index[0]
                elif sort_type == 'highest overall VAF':
                    max_vaf_sample = cluster_data.sum().sort_values(ascending=False).index[0]

                cluster_indexes = cluster_data.index.tolist()
                max_vaf_sample_index = cluster_data.columns.get_loc(max_vaf_sample)
                vaf_array = cluster_data.values.tolist()
                # streamlit.write(vaf_array)
                descending_vaf_array = sorted(vaf_array, key=lambda x: x[max_vaf_sample_index], reverse=True)
                # streamlit.write(descending_vaf_array)
                peak_sorted_vaf_array = []
                prepend = False
                for row in descending_vaf_array:
                    if prepend:
                        peak_sorted_vaf_array = [row] + peak_sorted_vaf_array
                        prepend = False
                    else:
                        peak_sorted_vaf_array = peak_sorted_vaf_array + [row]
                        prepend = True
                    # streamlit.write(peak_sorted_vaf_array)
                # streamlit.write(peak_sorted_vaf_array)
                cluster_data = pandas.DataFrame(peak_sorted_vaf_array, columns=vaf_columns, index=cluster_indexes)
                # streamlit.write(peak_sorted_vaf_array)

            before_df = pandas.DataFrame(index=list(range(data_frame.first_valid_index(), cluster_data.first_valid_index())), columns=vaf_columns)
            after_df = pandas.DataFrame(index=list(range(cluster_data.last_valid_index() + 1, data_frame.last_valid_index() + 1)), columns=vaf_columns)
            cluster_dataframes[cluster_colors[i]] = pandas.concat([before_df, cluster_data, after_df], axis=0)

        figure = go.Figure(data=[
            go.Surface(
                z=cluster_dataframes[color],
                x=sample_names,
                y=list(range(len(cluster_dataframes[color]))),
                colorscale=[[0, '#d9d9d9'], [1, color]],
                showscale=False
            ) for color in cluster_dataframes
        ])

        # To keep layout (rotate, zoom, etc. between plot edits) - see: https://stackoverflow.com/questions/68798315/how-to-update-plotly-plot-and-keep-ui-settings
        figure['layout']['uirevision'] = 'some string'
        
        figure.update_layout(
            width=1500,
            height=1500,
            scene=dict(
                yaxis=dict(
                    title='Variants',
                    showticklabels=False
                    # ticktext=list(range(len(data_frame))),
                    # nticks=len(data_frame)
                ),
                xaxis=dict(
                    title='Sample'
                ),
                zaxis=dict(
                    title='VAF'
                )
            )
        )

        streamlit.plotly_chart(figure, theme='streamlit', use_container_width=False)
    elif plot_type == 'ctDNA plot':
        from streamlit_sortables import sort_items

        # plot_theme = streamlit.sidebar.selectbox('Plot theme', ['plotly', 'plotly_white', 'plotly_dark', 'ggplot2', 'seaborn', 'simple_white', 'presentation', 'xgridoff', 'ygridoff', 'gridon', 'none'])
        # frame_width = streamlit.sidebar.slider('Frame width', min_value=1000, max_value=5000, step=100, value=3000)
        # frame_height = streamlit.sidebar.slider('Frame height', min_value=500, max_value=2500, step=100, value=900)

        # User input - target sample
        plasma_sample = streamlit.sidebar.selectbox('Plasma sample', [None] + sample_names)

        if not plasma_sample:
            streamlit.write('Please choose plasma, primary tumour, and metastasis samples in the menu to the left to see the plot.')

        else:

            comparison_samples_list = sample_names[:]
            comparison_samples_list.remove(plasma_sample)
            primary_tumor_samples_expander = streamlit.sidebar.expander('Primary tumor samples selector')
            metastatic_samples_expander = streamlit.sidebar.expander('Metastasis samples selector')
            primary_tumor_samples_dict = {}
            metastasis_samples_dict = {}

            for sample in comparison_samples_list:
                primary_tumor_samples_dict[sample] = primary_tumor_samples_expander.checkbox(sample, value=False, key=f'primary_option_{sample}')
                metastasis_samples_dict[sample] = metastatic_samples_expander.checkbox(sample, value=False, key=f'metastasis_option_{sample}')
            primary_tumor_samples = []
            metastasis_samples = []
            for sample in primary_tumor_samples_dict:
                if primary_tumor_samples_dict[sample]:
                    primary_tumor_samples.append(sample)
            for sample in metastasis_samples_dict:
                if metastasis_samples_dict[sample]:
                    metastasis_samples.append(sample)
            if primary_tumor_samples and metastasis_samples:

                # In the following, we want to assess VAF in primary tumor and metastases of mutations found in the plasma sample

                data = {
                    'VAF': [],
                    'ctDNA mutations': []
                }

                # User input - Include only ctDNA mutations
                ct_dna_mutations_only = streamlit.sidebar.checkbox('Include only ctDNA mutations', value=False)
                ct_dna_data_frame = data_frame[data_frame[f'VAF_{plasma_sample}'] > 0]
                if ct_dna_mutations_only:
                    test_df = ct_dna_data_frame
                else:
                    test_df = data_frame

                shared_n = 0
                primary_n = 0
                metastasis_n = 0
                for index, row in test_df.iterrows():
                    found_in_primary = False
                    found_in_metastasis = False
                    for sample in primary_tumor_samples:
                        if row[f'VAF_{sample}'] > 0:
                            found_in_primary = True
                            break
                    for sample in metastasis_samples:
                        if row[f'VAF_{sample}'] > 0:
                            found_in_metastasis = True
                            break
                    if found_in_primary and found_in_metastasis:
                        data['VAF'].append(row[f'VAF_{plasma_sample}'])
                        data['ctDNA mutations'].append('Shared')
                        shared_n += 1
                    elif found_in_primary:
                        data['VAF'].append(row[f'VAF_{plasma_sample}'])
                        data['ctDNA mutations'].append('Primary-specific')
                        primary_n += 1
                    elif found_in_metastasis:
                        data['VAF'].append(row[f'VAF_{plasma_sample}'])
                        data['ctDNA mutations'].append('Metastasis-specific')
                        metastasis_n += 1

                plotting_data = pandas.DataFrame.from_dict(data)
                # streamlit.write(plotting_data)

                # Define Expander for visual appearance
                visual_appearance = streamlit.sidebar.expander('Edit visual appearance')

                plot_theme = visual_appearance.selectbox('Plot theme', ['plotly', 'plotly_white', 'plotly_dark', 'ggplot2', 'seaborn', 'simple_white', 'presentation', 'xgridoff', 'ygridoff', 'gridon', 'none'])

                frame_layout = streamlit.sidebar.expander('Frame layout')
                frame_layout.write('If any of the plots are "cut off", you can increase the size of the plot display frame here - rarely necessary.')
                frame_width = frame_layout.slider('Frame width', min_value=1000, max_value=5000, step=100, value=3000)
                frame_height = frame_layout.slider('Frame height', min_value=500, max_value=2500, step=100, value=900)

                # User input - X-axis sorting
                category_order = ['Primary-specific', 'Metastasis-specific', 'Shared']
                order_widget = visual_appearance.container()
                order_widget.write('Choose x-axis order by dragging:')
                with order_widget:
                    sorted_x_items = sort_items(category_order)

                # User input - plot type
                plot_type = visual_appearance.radio('Plot type', ['Box', 'Violin', 'Combined'])

                # User input - plot size
                plot_height = visual_appearance.number_input('Plot height', min_value=200, max_value=1000, value=600, step=20)
                plot_width = visual_appearance.number_input('Plot width', min_value=100, max_value=1000, value=400, step=20)

                # User input - edit y-axis
                edit_y_axis = visual_appearance.checkbox('Edit Y-axis range', value=False)

                # User input - y-axis range
                yaxis_max = None
                yaxis_min = None
                if edit_y_axis:
                    visual_appearance.write('Please define values for both min and max to see an effect.')
                    yaxis_max = visual_appearance.number_input('Y-axis maximum', min_value=0.1, max_value=1.0, value=None)
                    yaxis_min = visual_appearance.number_input('Y-axis minimum', min_value=-1.0, max_value=0.00, value=0.0)

                if plot_type == 'Box':
                    figure = px.box(plotting_data, y='VAF', color='ctDNA mutations', width=plot_width, height=plot_height, template=plot_theme, category_orders={'ctDNA mutations': sorted_x_items})
                elif plot_type == 'Violin':
                    # Make violin plot from the data
                    figure = px.violin(plotting_data, y='VAF', color='ctDNA mutations', width=plot_width, height=plot_height, template=plot_theme, category_orders={'ctDNA mutations': sorted_x_items})
                elif plot_type == 'Combined':
                    figure = px.violin(plotting_data, y='VAF', color='ctDNA mutations', width=plot_width, height=plot_height, box=True, template=plot_theme, category_orders={'ctDNA mutations': sorted_x_items})

                if yaxis_max is not None and yaxis_min is not None:
                    figure.update_layout(
                        yaxis_range=[yaxis_min, yaxis_max],
                    )

                # Horizontal legends on op
                # figure.update_layout(legend=dict(
                #     orientation="h",
                #     yanchor="bottom",
                #     y=1.02,
                #     xanchor="right",
                #     x=1
                # ))
                # figure.update_layout(
                #     hoverlabel_align='left'  # Necessary for streamlit to make text for all labels align left
                # )

                # streamlit.text(f'ctDNA mutations ratio:\t\t\t{round(len(ct_dna_data_frame) / len(data_frame), 3)}')
                streamlit.text(f'# Primary-specific mutations:\t\t{primary_n}')
                streamlit.text(f'# Metastasis-specific mutations:\t{metastasis_n}')
                streamlit.text(f'# Shared mutations:\t\t\t{shared_n}')

                # streamlit.plotly_chart(figure, theme='streamlit', use_container_width=False)  # Old way
                streamlit_components_v1.html(figure.to_html(), width=frame_width, height=frame_height, scrolling=True)

                export_to_svg_ui(figure)

    elif plot_type == 'Fish plot':
        class Node:

            def __init__(self, node_id, parent_id=None, child_id=None):
                self.node_id = node_id
                self.parent_id = parent_id
                self.level = None
                self.children_ids = {child_id} if child_id is not None else set()

            def add_parent(self, parent_id):
                self.parent_id = parent_id

            def add_child(self, child_id):
                self.children_ids.add(child_id)

            def set_level(self, level):
                self.level = level

            def has_children(self):
                return len(self.children_ids) > 0


        class Tree:
            """
            The purpose of this class is to get depth and level for each node and for the tree
            as this is needed for plotting. It depends on the tree branches being connected, but whether this is the case or not is
            not controlled. As it is designed to process CONIPHER generated trees, it is assumed to hold true.
            """

            def __init__(self):
                self.nodes = {}  # Node ID > node object dict
                self.no_children_nodes_ids = set()
                self.parent_nodes_ids = []
                self.children_nodes_ids = []
                self.top_node_id = None
                self.levels = {}  # Contains lists of node IDs in each level
                self.height = 0  # Trunk is 1-indexed

            def add_levels(self, node_id, level):
                """
                Adds given level to node with given id and does the same recursively for each child while incrementing the
                level one for each recursion.

                :param level:
                :param node_id:
                :return:
                """
                node = self.nodes[node_id]
                node.set_level(level)
                self.levels.setdefault(level, []).append(node_id)
                for child_node_id in node.children_ids:
                    self.add_levels(child_node_id, level + 1)


            def order_nodes_according_to_tree(self, parent_node_id, descendants):
                """
                The engine of get_tree_order_of_nodes. Runs through the tree recursively and adds node IDs to
                self.node_order (defined in get_tree_order_of_nodes) if they are in descendants
                :param parent_node_id:
                :param descendants:
                :return:
                """
                if parent_node_id in descendants:
                    self.node_order.append(parent_node_id)
                parent_node = self.nodes[parent_node_id]
                for child_node_id in parent_node.children_ids:
                    self.order_nodes_according_to_tree(child_node_id, descendants)


            def get_tree_order_of_nodes(self, parent_node_id, descendants):
                """
                Given a parent node and a list of IDs of some or all descendants of that node, this function returns the
                IDs of the parent node and the given descendants in order relative to the tree structure starting from
                the top node as opposed to numerical order.
                :param parent_node_id:
                :param descendants:
                :return: List of parent and descendant node IDs in tree order
                """
                self.node_order = []
                self.order_nodes_according_to_tree(parent_node_id, descendants)
                return self.node_order


            def add_branches(self, branches):
                """
                Adds given branches to the tree and calculates top node, tree height and adds level information to the nodes
                :param branches:    A list of tab delimited two-node-branches - corresponding to CONIPHER output
                :return:
                """
                for branch in branches:
                    branch_nodes = branch.split('-')
                    parent_node_id = branch_nodes[0]
                    child_node_id = branch_nodes[1]
                    self.parent_nodes_ids.append(parent_node_id)
                    self.children_nodes_ids.append(child_node_id)
                    if not parent_node_id in self.nodes:

                        self.nodes[parent_node_id] = Node(parent_node_id, child_id=child_node_id)
                    else:
                        self.nodes[parent_node_id].add_child(child_node_id)
                        if parent_node_id in self.no_children_nodes_ids:
                            self.no_children_nodes_ids.remove(parent_node_id)
                    if not child_node_id in self.nodes:

                        self.nodes[child_node_id] = Node(child_node_id, parent_id=parent_node_id)
                        self.no_children_nodes_ids.add(child_node_id)
                    else:
                        self.nodes[child_node_id].add_parent(parent_node_id)


                # Get top node
                self.top_node_id = [x for x in self.parent_nodes_ids if x not in self.children_nodes_ids][0]

                # Add levels starting from the top node
                self.add_levels(self.top_node_id, 1)

                # Calculate tree height
                end_node_levels = []
                for end_node_id in self.no_children_nodes_ids:
                    end_node_levels.append(self.nodes[end_node_id].level)

                self.height = max(end_node_levels)


        number_of_plots = len(clone_proportions)
        samples = list(clone_proportions.keys())

        if number_of_plots > 1:

            # User input - Which samples to plot
            sample = streamlit.sidebar.selectbox('Please select which samples to plot', samples)
        else:
            sample = samples[0]

        tree = Tree()
        tree.add_branches(tree_branches)

        # Construct pyfish tables using the Tree object
        branches = {'ParentId': tree.parent_nodes_ids, 'ChildId': tree.children_nodes_ids}
        parent_tree_df = pandas.DataFrame.from_dict(branches)

        sample_proportion = clone_proportions[sample]

        represented_clones = {}  # Contains clones that have a proportion > 0
        if tree.height >= 2:
            population = {
                'Id': [],
                'Step': [],
                'Pop': [],
                'cluster_id': []
            }
            for level in list(range(1, tree.height + 1)):
                for node_id in tree.levels[level]:
                    color_index = int(node_id) - 1
                    if sample_proportion[node_id] == 'NA':
                        proportion = 0.0
                    else:
                        proportion = float(sample_proportion[node_id])

                    if proportion > 0.0:
                        represented_clones[node_id] = {
                            'display_id': f'clone_{node_id}',
                            'proportion': proportion,
                            'color': cluster_colors[color_index]
                        }

                    # We define just two steps - the starting point (step == 0) and the end point (step == tree height)

                    # For all clones,we assume they have zero representation at step 0
                    population['Id'].append(int(node_id))
                    population['Step'].append(level - 1)
                    population['Pop'].append(0.0)
                    population['cluster_id'].append(color_index)

                    population['Id'].append(int(node_id))
                    population['Step'].append(tree.height)
                    population['Pop'].append(proportion)
                    population['cluster_id'].append(color_index)

        else:
            pass

        populations_df = pandas.DataFrame.from_dict(population)

        parent_tree_df = parent_tree_df.astype(int)



        # The documentation in MatplotLib and pyfish for making and using custom color maps is not very good. According to
        # the mpl documentation, ListedColormap requires a list of colors as either color names or RGB(A) values. However,
        # RGBA values as tuples seem to fail (perhaps they should be written as a string). Hex values seem to work, though.
        # https://matplotlib.org/stable/api/_as_gen/matplotlib.colors.ListedColormap.html#matplotlib.colors.ListedColormap
        #
        # In the call to process_data below (PyFish), the named argument, cmap_name, seems to imply a string. However, it
        # must be the color map object returned by ListedColormap. Even, using the name argument in the call to
        # ListedColormap and using that name does not work.
        #
        # Lastly, the colormap returned by ListedColormap seems to be one-indexed, so cluster-IDs, which are also
        # one-indexed can be used directly.

        # User input - figure width and height
        # Define Expander for visual appearance
        visual_appearance = streamlit.sidebar.expander('Edit Fish plot')

        # User input - figure width
        figure_width = visual_appearance.slider('Figure width', min_value=300, max_value=3000, step=10, value=1180)

        # User input - figure width
        figure_height = visual_appearance.slider('Figure height', min_value=150, max_value=2000, step=10, value=800)

        # User input - interpolation
        smooth = visual_appearance.slider('Smoothing', min_value=0, max_value=10, step=1, value=6)

        # User input - Proportions position
        proportions_horizontal = visual_appearance.slider('Proportions - horizontal position', min_value=5, max_value=100, step=1, value=32)
        proportions_vertical = visual_appearance.slider('Proportions - vertical position', min_value=0, max_value=1000, step=1, value=100)

        # User input proportions sorting
        sort_proportions = visual_appearance.radio('Sort proportions by:', ['Nummerical clone ID', 'Phylogenetic tree order'])

        # Because PyFish maps colors relative to max and min values of the cluster_color column, we must make sure that
        # the used list of colors match the range of those. Otherwise, the coloring is off in inconsistent manners.
        # We therefore find numerical min and max clones and use those to slice the cluster_colors list
        min_clone = min(parent_tree_df['ParentId'].min(), parent_tree_df['ChildId'].min()) - 1
        max_clone = max(parent_tree_df['ParentId'].max(), parent_tree_df['ChildId'].max()) - 1
        my_matplotlib_color_map = ListedColormap(cluster_colors[min_clone: max_clone + 1])
        data = process_data(populations_df, parent_tree_df, interpolation=1, smooth=smooth, absolute=True, cmap_name=my_matplotlib_color_map, color_by='cluster_id')
        setup_figure(width=figure_width, height=figure_height, absolute=True)
        fish_plot(*data)

        fig = matplotlib_pyplot.gcf()

        # Does not work - mpl_to_plotly was experimental and is no longer supported
        # figure = plotly_tools.mpl_to_plotly(fig)
        #
        # streamlit_components_v1.html(figure.to_html(), width=figure_width, height=figure_height, scrolling=True)
        # export_to_svg_ui(figure)

        html_fig = mpld3.fig_to_html(fig)

        cluster_proportions_html = f"""
        <style>
            #cluster_proportions {{
                position: relative;
                top: {proportions_vertical}px;
            }}
          .clone_color {{
            width: 12px;
            height: 12px;
            display: block;
            float: left;
            border-radius: 20px;
            margin-right: 10px;
            position: relative;
            top: 6px;
          }}
        </style>
        <div id="cluster_proportions">
            <h5>
            Proportions
            </h5>

        """
        if sort_proportions == 'Nummerical clone ID':
            represented_clones = dict(sorted(represented_clones.items(), key=lambda clone: int(clone[0])))
        elif sort_proportions == 'Phylogenetic tree order':
            highest_represented_clone_id = sorted(list(represented_clones.keys()), key=lambda clone_id: tree.nodes[clone_id].level)[0]
            represented_clone_ids = set(represented_clones.keys())
            represented_clones_in_tree_order = tree.get_tree_order_of_nodes(highest_represented_clone_id, represented_clone_ids)
            mapping = {elm: i for i, elm in enumerate(represented_clones_in_tree_order)}
        for clone in represented_clones:
            cluster_proportions_html += f"""
            <p>
                <span class="clone_color" style="background-color: {represented_clones[clone]['color']}"></span>
                <span>{represented_clones[clone]['display_id']} ({represented_clones[clone]['proportion']} %)</span>
            </p>
            """
        cluster_proportions_html += '</div>'

        # Separate plot and Clone proportions display in two columns
        col1, col2, col3, col4 = streamlit.columns([0.4, 0.1, 0.2, 0.3])
        horizontal = proportions_horizontal / 100
        col1, col2 = streamlit.columns([horizontal, 1 - horizontal])
        with col1:
            streamlit_components_v1.html(html_fig, width=figure_width, height=figure_height, scrolling=False)
        with col2:
            streamlit.html(cluster_proportions_html)


