# -*- coding: utf-8 -*-
import pandas as pd

import dash
#from dash import Dash
from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html
import colorlover as cl
import numpy as np

df_targets = pd.read_csv('data_out.csv')


def create_traces(df_targets, highlight_target, min_val, max_val):
    traces = []
    bupu = cl.scales['9']['div']['PiYG']
    scale = cl.interp(bupu, 500)  # Map color scale to 500 bins
    grouped = df_targets.groupby('target_chemblid')
    for chembl_id, df_by_target in grouped:
        # total = df_by_target['cumsum'].max()
        # if total <50:
        #     continue
        target = df_by_target['target'].iloc[0]
        growth = df_by_target['best_phase'].max()
        peak_to_current = df_by_target['peak_to_current'].max()
        # if peak_to_current <=min_val or peak_to_current >=max_val:
        #     continue
        color = scale[int(peak_to_current * 500 - 1)]
        x = df_by_target['year']

        y = df_by_target['r_difference']

        if growth >= 4:
            # color = scale[0]
            dash = None
            legendgroup = '4'
            name = 'Phase 4'
        elif growth >= 3:
            # color = scale[1]
            dash = 'dot'
            legendgroup = '3'
            name = 'Phase 3'
        elif growth >= 2:
            # color = scale[2]
            legendgroup = '2'
            dash = 'dashdot'
            name = 'Phase 2'
        elif growth >= 1:
            # color = scale[3]
            legendgroup = '1'
            dash = None
            name = 'Phase 1'
        else:
            # color = scale[4]
            legendgroup = '0'
            name = 'Discovery'
            dash = 'dash'

        if highlight_target and target not in highlight_target:
            color = 'lightgrey'

        hoverinfo = 'text'
        width = 1
        if highlight_target and target in highlight_target:
            width = 2.5
            if color == 'lightgrey':
                color = 'grey'
            hoverinfo = 'text'
        elif highlight_target:
            width = 0.5
            hoverinfo = 'none'

        traces.append({
            'x': x,
            'y': y,
            'mode': 'lines',
            'line': {
                'color': color,
                'width': width,
                'dash': None,
                'shape': 'spline'
            },
            'text': ["{},{},{},{}".format(target, chembl_id, b, c) for b, c in zip(x, y)],
            'legendgroup': legendgroup,
            'name': name,
            'customdata': (target,chembl_id),
            'hoverinfo': hoverinfo,
            'showlegend': (
                False if legendgroup in [t['legendgroup'] for t in traces]
                else True
            )
        })
    return traces


def create_figure(df_targets, highlight_target=None, skip_labels=[], show_only=[], min_val=0, max_val=1):
    # min_year = df_targets['YEAR'].min()
    # max_year = df_targets['YEAR'].max()
    # Construct the charts
    # scale = cl.flipper()['seq']['5']['Reds']

    traces = create_traces(df_targets, highlight_target, min_val, max_val)
    annotations = []
    # reorder traces to reorder legend items



    if not highlight_target:
        def get_trace_index(traces, legendgroup):
            for i, trace in enumerate(traces):
                if trace['showlegend'] and trace['legendgroup'] == legendgroup:
                    return i

        try:
            traces.insert(0, traces.pop(get_trace_index(traces, '0')))
            traces.insert(0, traces.pop(get_trace_index(traces, '1')))
            traces.insert(0, traces.pop(get_trace_index(traces, '2')))
            traces.insert(0, traces.pop(get_trace_index(traces, '3')))
            traces.insert(0, traces.pop(get_trace_index(traces, '4')))
        except:
            pass
    else:
        # move highlighted traces to the end
        for i, trace in enumerate(traces):
            if trace['line']['width'] != 2.5:
                traces.insert(0, traces.pop(i))

    if not highlight_target:
        annotations = [{
            'x': 0.8, 'xref': 'paper', 'xanchor': 'left',
            'y': 0.95, 'yref': 'paper', 'yanchor': 'bottom',
            'text': '<b>Max Phase Achieved</b>',
            'showarrow': False
        }]

    layout = {
        'xaxis': {
            'showgrid': False,
            'showline': False,
            'zeroline': False,
            'showticklabels': False,
            'ticks': '',
            'title': 'Year'
        },
        'yaxis': {
            'showgrid': False,
            'showticklabels': False,
            'zeroline': False,
            # 'type':'log',
            'ticks': '',
            'title': 'Relative number of compounds published per year'
        },
        'showlegend': True,  # not bool(highlight_target),
        'hovermode': 'closest',
        'legend': {
            'x': 0.8,
            'y': 0.95,
            'xanchor': 'left'
        },
        'annotations': annotations,
        'margin': {'t': 20, 'b': 20, 'r': 0, 'l': 20},
        'font': {'size': 12}
    }

    return {'data': traces, 'layout': layout}


def serve_layout():
    layout = html.Div([

        dcc.Markdown('''
        # Target Attention Over Time
    
        The number of compounds published for a target in a given year gives an indication of the confidence at that 
        point in time that the target may yield a successful drug. The following plots show the moving average of the 
        relative number of compounds published each year since 1990 for targets in ChEMBL. Each line is coloured by the
        target's *recent attention score*, the ratio between its current and peak publication rate.
        
        ***
    
            '''.replace('  ', ''), className='container',
                     containerProps={'style': {'maxWidth': '650px'}}),

        dcc.Markdown('''
        ## Contemporary Projects

       These targets represent those for which compounds are being discovered at or near their peak rate since 1990. 
       They can be filtered by the maximum phase achieved for the target by double clicking on the desired phase in the 
       legend. Filtering to show only phase 4 suggests that majority of contemporary projects focus on known drug targets. 

                '''.replace('  ', ''), className='container',
                     containerProps={'style': {'maxWidth': '650px'}}),

        html.Div([
            html.Div([
                dcc.Graph(
                    figure=create_figure(df_targets[
                                             (df_targets['peak_to_current'] > 0.8) & (
                                             df_targets['peak_to_current'] < 2)]),
                    id='contemporary',
                    style={'height': '85vh'})
            ], className='eight columns'),
            html.Div([
                html.Iframe(
                    # enable all sandbox features
                    # see https://developer.mozilla.org/en-US/docs/Web/HTML/Element/iframe
                    # this prevents javascript from running inside the iframe
                    # and other things security reasons
                    id='chembl_widget_contemporary',
                    src="http://chembl-glados.herokuapp.com/target_report_card/CHEMBL204/embed/approved_drugs_clinical_candidates/",
                    style={'height': '80vh', 'width':'100%', 'border':'none'}
                )
            ], className='four columns'),
        ], className='row',  style={'padding-bottom':'10vh', 'padding-top':'10vh'}),

        dcc.Markdown('''
        ## Former Projects

        Compounds are currently being discovered at a much lower rate for these targets than in their peak years.
        This could indicate that a suitable drug has been found, and it is no longer beneficial to persue the target, 
        or that it has since been found that it is not possible to find a drug for the target.  
        In order to differentiate between these two scenarios, the targets can be filtered by their maximum phase
        by clicking on the legend. For example, double clicking 'Discovery' will only display targets that have
        yet to make it into the clinic. These targets represent a potential set of intractable targets. 

            '''.replace('  ', ''), className='container',
                     containerProps={'style': {'maxWidth': '650px'}}),

        html.Div([
            html.Div([
                html.Iframe(
                    # enable all sandbox features
                    # see https://developer.mozilla.org/en-US/docs/Web/HTML/Element/iframe
                    # this prevents javascript from running inside the iframe
                    # and other things security reasons
                    id='chembl_widget_former',
                    src="http://chembl-glados.herokuapp.com/target_report_card/CHEMBL204/embed/approved_drugs_clinical_candidates/",
                    style={'height': '80vh', 'width':'100%', 'border':'none'}
                )
            ], className='four columns'),

            html.Div([
                dcc.Graph(
                    figure=create_figure(df_targets[
                                             (df_targets['peak_to_current'] > 0) & (
                                             df_targets['peak_to_current'] < 0.3)]),
                    id='former',
                    style={'height': '85vh'})
            ], className='eight columns'),
        ],className='row',  style={'padding-bottom':'10vh','padding-top':'10vh'}),

        dcc.Markdown('''
    
        ## Interactive Overview 
        
        Targets with green lines show a high ratio, meaning that compounds are bing published at a similar rate
        to the peak years. Targets with pink lines show a low ratio, meaning that far fewer compounds have been 
        published in recent years compared to the target's peak years. Individual targets can be selected by 
        clicking on the graph, or selecting in the dropdown box. Custom ranges for the Peak:Recent values can be
        selected using the range slider.
        *** 
      
        '''.replace('  ', ''), className='container',
                     containerProps={'style': {'maxWidth': '650px'}}),

        html.Div([
            html.Div([
                dcc.Dropdown(
                    options=[
                        {'label': c, 'value': c}
                        for c in sorted(list(df_targets.target.unique()))
                    ],
                    value=[None, None],
                    multi=True,
                    id='category-filter'
                )
            ], className='row', style={'padding-bottom': '10px'}),
            html.Div([
                html.Label('Current attention score', id='num_compounds_label', className='inline'),
                dcc.RangeSlider(
                    id='peak_current',
                    min=0,
                    max=1,
                    step=0.01,
                    marks={round(x * 0.1, 1): round(x * 0.1, 1) for x in range(0, 11)},
                    value=[0, 1]

                )], className='eight columns'),

            html.Div([
                html.Button(id='subbutton', n_clicks=0, children='Update Graph')
            ], className='three columns'),
        ], className='container', style={'maxWidth': '650px'}),
        html.Div([
            dcc.Graph(
                figure=create_figure(df_targets), id='overview',
                style={'height': '85vh'})
        ]),
    ])

    return layout


#server = Flask(__name__)
#app = Dash(__name__, server=server, csrf_protect=False)
app = dash.Dash(__name__)
server = app.server


app.css.append_css({
    'external_url': (
        'https://cdn.rawgit.com/chriddyp/0247653a7c52feb4c48437e1c1837f75'
        '/raw/a68333b876edaf62df2efa7bac0e9b3613258851/dash.css'
    )

})

app.css.append_css({
    'external_url': (
        "https://codepen.io/mikesmith1611/pen/QOKgpG.css"
    )

})

app.layout = serve_layout()


@app.callback(
    Output('overview', 'figure'),
    [Input('subbutton', 'n_clicks'), Input('overview', 'clickData')],
    [State('peak_current', 'value'), State('category-filter', 'value')])
def filter(n_clicks, selected_value, range_value, category_filter):
    if selected_value is None:
        filter_li = category_filter
        filter_li = [f for f in filter_li if f is not None]
        df_filtered = df_targets[
            (df_targets['peak_to_current'] > range_value[0]) & (df_targets['peak_to_current'] < range_value[1])]
    else:
        target = selected_value['points'][0]['text'].split(',')[0]
        filter_li = [target]
        df_filtered = df_targets
    figure = create_figure(
        df_filtered,
        filter_li if len(filter_li) > 0 else None,
        skip_labels=['-'],
        min_val=range_value[0],
        max_val=range_value[1]
    )

    for trace in figure['data']:
        trace['hoverinfo'] = 'text'

    return figure

@app.callback(
    Output('chembl_widget_contemporary', 'src'),
    [Input('contemporary', 'clickData')],
    )
def get_contemporary_widget(selected_value):

    print selected_value
    chembl_id = selected_value['points'][0]['text'].split(',')[1]
    src = "http://chembl-glados.herokuapp.com/target_report_card/{}/embed/approved_drugs_clinical_candidates/".format(chembl_id)

    return src

@app.callback(
    Output('chembl_widget_former', 'src'),
    [Input('former', 'clickData')],
    )
def get_former_widget(selected_value):

    print selected_value
    chembl_id = selected_value['points'][0]['text'].split(',')[1]
    src = "http://chembl-glados.herokuapp.com/target_report_card/{}/embed/approved_drugs_clinical_candidates/".format(chembl_id)

    return src

if __name__ == '__main__':
    app.run_server(port=9999)
