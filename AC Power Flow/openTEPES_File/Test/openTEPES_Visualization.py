# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.33 - June 15, 2020

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import networkx as nx
import plotly.graph_objs as go
from decimal import Decimal

import pandas as pd
from colour import Color
from datetime import datetime, timedelta, date
from textwrap import dedent as d
import json
import re

# import the css template, and pass the css template into dash
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app.title = "SNetwork"

CaseName  = 'RTS96'

InputNode = "N_201"
InputDate = "2030-03-10 01:00:00"
InputTime = 18

##############################################################################################################################################################
def network_graph(InputDate, InputTime):
#     InputDate                  = datetime.strptime(re.split('T| ', InputDate)[0], '%Y-%m-%d %H:%M:%S')
    InputDate                  = datetime.strptime(re.split('T| ', InputDate)[0], '%Y-%m-%d')
    InputDate                  = datetime(InputDate.year, InputDate.month, InputDate.day, 0, 0, 0) + timedelta(hours=InputTime)
    

    dfNetwork                  = pd.read_csv('oT_Result_NetworkUtilization_ToPT_' + CaseName + '.csv')
    dfNetworkFlow              = pd.read_csv('oT_Result_NetworkFlow_ToPT_' + CaseName + '.csv')
    dfDemand                   = pd.read_csv('oT_Data_Demand_' + CaseName + '.csv', index_col=[0, 1, 2])
    dfNodeLocation             = pd.read_csv('oT_Data_NodeLocation_' + CaseName + '.csv', index_col=[0])
    dfGeneration               = pd.read_csv('oT_GenDict.csv', index_col=[0    ])
    dfNode                     = pd.read_csv('oT_NodeDict.csv', index_col=[0    ])
    dfGenOut                   = pd.read_csv('oT_Result_GenerationOutput_ToPT_RTS96.csv', index_col=[0,1,2,3,4])
    dfPNS                      = pd.read_csv('oT_Result_NetworkPNS_ToPT_RTS96.csv')

    # DataFrame stack
    dfDemand                   = dfDemand.stack()
    # To set index
    dfDemand.index.names       = ['Scenario', 'Period', 'LoadLevel', 'Node']
    dfGenOut.index.names       = ['Scenario', 'Period', 'LoadLevel', 'GenUnit', 'MW']
    # Reset Indexes
    dfDemand                   = dfDemand.reset_index()
    dfGenOut                   = dfGenOut.reset_index()
    dfGenOut                   = dfGenOut.assign(Technology = dfGenOut['GenUnit'])
    dfGenOut                   = dfGenOut.assign(Node       = dfGenOut['GenUnit'])
    dfGenOut                   = dfGenOut.assign(System     = dfGenOut['GenUnit'])

    for i in dfGeneration.index:
        dfGenOut.loc[dfGenOut['Technology']   == i, 'Technology']   = dfGeneration['Technology'][i]
        dfGenOut.loc[dfGenOut['Node']         == i, 'Node']         = dfGeneration['Node'][i]
        dfGenOut.loc[dfGenOut['System']       == i, 'System']       = dfGeneration['System'][i]
    # Rename columns
    dfNetwork                  = dfNetwork.rename(columns={"Scenario": "Scenario", "Period": "Period"
                                                           , "LoadLevel": "LoadLevel", "InitialNode": "InitialNode"
                                                           , "FinalNode": "FinalNode", "Circuit": "Circuit", "0": "MW"})
    dfNetworkFlow              = dfNetworkFlow.rename(columns={"Scenario": "Scenario", "Period": "Period"
                                                           , "LoadLevel": "LoadLevel", "InitialNode": "InitialNode"
                                                           , "FinalNode": "FinalNode", "Circuit": "Circuit", "vFlow": "MW"})
    dfDemand                   = dfDemand.rename(columns={"Scenario": "Scenario", "Period": "Period"
                                                          , "LoadLevel": "LoadLevel", "Node": "Node", 0: "MW"})
    
    # Rename columns
    dfPNS                      = dfPNS.rename(columns={"Scenario": "Scenario", "Period": "Period", "LoadLevel": "LoadLevel"
                                                       , "Node": "GenUnit", "vENS": "MW"})
    dfPNS                      = dfPNS.assign(Technology = dfPNS['GenUnit'])
    dfPNS                      = dfPNS.assign(Node       = dfPNS['GenUnit'])
    dfPNS                      = dfPNS.assign(System     = dfPNS['GenUnit'])
    dfPNS                      = dfPNS.assign(GraphType  = dfPNS['GenUnit'])

    for i in dfNode.index:
        dfPNS.loc[dfPNS['Technology']   == i, 'Technology']   = dfNode['Technology'][i]
        dfPNS.loc[dfPNS['Node']         == i, 'Node']         = dfNode['Node'][i]
        dfPNS.loc[dfPNS['System']       == i, 'System']       = dfNode['System'][i]
        dfPNS.loc[dfPNS['GraphType']    == i, 'GraphType']    = dfNode['GraphType'][i]
    # Change data format
    dfNetwork['LoadLevel']     = pd.to_datetime(dfNetwork.LoadLevel)
    dfNetworkFlow['LoadLevel'] = pd.to_datetime(dfNetworkFlow.LoadLevel)
    dfDemand['LoadLevel']      = pd.to_datetime(dfDemand.LoadLevel)
    dfGenOut['LoadLevel']      = pd.to_datetime(dfGenOut.LoadLevel)
    dfPNS['LoadLevel']         = pd.to_datetime(dfPNS.LoadLevel)
    # Add Columns datetime
    dfNetwork                  = dfNetwork.assign(Year=dfNetwork['LoadLevel'].dt.year)
    dfNetwork                  = dfNetwork.assign(Month=dfNetwork['LoadLevel'].dt.month)
    dfNetwork                  = dfNetwork.assign(Day=dfNetwork['LoadLevel'].dt.day)
    dfNetwork                  = dfNetwork.assign(Hour=dfNetwork['LoadLevel'].dt.hour)
                 
    dfNetworkFlow              = dfNetworkFlow.assign(Year=dfNetworkFlow['LoadLevel'].dt.year)
    dfNetworkFlow              = dfNetworkFlow.assign(Month=dfNetworkFlow['LoadLevel'].dt.month)
    dfNetworkFlow              = dfNetworkFlow.assign(Day=dfNetworkFlow['LoadLevel'].dt.day)
    dfNetworkFlow              = dfNetworkFlow.assign(Hour=dfNetworkFlow['LoadLevel'].dt.hour)
             
    dfDemand                   = dfDemand.assign(Year=dfDemand['LoadLevel'].dt.year)
    dfDemand                   = dfDemand.assign(Month=dfDemand['LoadLevel'].dt.month)
    dfDemand                   = dfDemand.assign(Day=dfDemand['LoadLevel'].dt.day)
    dfDemand                   = dfDemand.assign(Hour=dfDemand['LoadLevel'].dt.hour)
                 
    dfGenOut                   = dfGenOut.assign(Year=dfGenOut['LoadLevel'].dt.year)
    dfGenOut                   = dfGenOut.assign(Month=dfGenOut['LoadLevel'].dt.month)
    dfGenOut                   = dfGenOut.assign(Day=dfGenOut['LoadLevel'].dt.day)
    dfGenOut                   = dfGenOut.assign(Hour=dfGenOut['LoadLevel'].dt.hour)
                 
    dfPNS                      = dfPNS.assign(Year=dfPNS['LoadLevel'].dt.year)
    dfPNS                      = dfPNS.assign(Month=dfPNS['LoadLevel'].dt.month)
    dfPNS                      = dfPNS.assign(Day=dfPNS['LoadLevel'].dt.day)
    dfPNS                      = dfPNS.assign(Hour=dfPNS['LoadLevel'].dt.hour)

    # filter the record by datetime, to enable interactive control through the input box
    dfNetwork                  = dfNetwork.loc[(dfNetwork['Year'] == InputDate.year) & (dfNetwork['Month'] == InputDate.month) 
                                               & (dfNetwork['Day'] == InputDate.day) & (dfNetwork['Hour'] == InputDate.hour)]
    
    dfNetworkFlow              = dfNetworkFlow.loc[(dfNetworkFlow['Year'] == InputDate.year) & (dfNetworkFlow['Month'] == InputDate.month) 
                                               & (dfNetworkFlow['Day'] == InputDate.day) & (dfNetworkFlow['Hour'] == InputDate.hour)]

    dfDemand                   = dfDemand.loc[(dfDemand['Year']    == InputDate.year) & (dfDemand['Month']  == InputDate.month) 
                                              & (dfDemand['Day']  == InputDate.day) & (dfDemand['Hour']  == InputDate.hour)]
    
    dfGenOut                   = dfGenOut.loc[(dfGenOut['Year']   == InputDate.year) & (dfGenOut['Month']  == InputDate.month) 
                                              & (dfGenOut['Day']  == InputDate.day) & (dfGenOut['Hour']  == InputDate.hour)]

    dfPNS                      = dfPNS.loc[(dfPNS['Year']         == InputDate.year) & (dfPNS['Month']     == InputDate.month) 
                                           & (dfPNS['Day']     == InputDate.day) & (dfPNS['Hour']     == InputDate.hour)]

    pDemand                    = dfDemand[['Node', 'MW']]
    pDemand                    = pDemand.set_index('Node')
    pPNS                       = dfPNS[['Node', 'MW']]
    pPNS                       = pPNS.set_index('Node')

    #Add Column Color Type
    dfNetwork                  = dfNetwork.assign(ColorType = dfNetwork['MW'])
    dfNetwork                  = dfNetwork.assign(vFlow=dfNetworkFlow['MW'])
    dfNetwork.to_csv('oT_Result_Output.csv', sep=',')
    for i in dfNetwork.index:
        if dfNetwork['ColorType'][i] > 0.75:
            dfNetwork['ColorType'][i] = 'lightcoral'
        else:
            dfNetwork['ColorType'][i] = 'darkgray'

    GenToNode                  = dfGenOut.groupby(['Node']).sum()['MW']
    
    dfGenToNode = pd.DataFrame(columns=['MW'])
    for i in range(len(dfNode.index)):
#     for i in 0:len(dfNode.index):
        dfGenToNode = dfGenToNode.append({'MW': 0}, ignore_index=True)

    dfGenToNode.index = list(dfNode.index)

    for i in GenToNode.index:
        dfGenToNode['MW'][i] = GenToNode[i]
    
    accountSet = set()
    for index in dfNetwork.index:
        accountSet.add(dfNetwork['InitialNode'][index])
        accountSet.add(dfNetwork['FinalNode'][index])

    # to define the centric point of the networkx layout
    shells=[]
    shell1=[]
    shell1.append(InputNode)
    shells.append(shell1)
    shell2=[]
    for ele in accountSet:
        if ele!=InputNode:
            shell2.append(ele)
    shells.append(shell2)


    G = nx.from_pandas_edgelist(dfNetwork, 'InitialNode', 'FinalNode', ['InitialNode', 'FinalNode', 'vFlow', 'MW', 'ColorType'], 
                                create_using=nx.MultiGraph())

    nx.set_node_attributes(G, dfGenToNode.to_dict(), 'GenOut')
    
    # We have to set the population attribute for each of the 14 nodes
    for i in list(G.nodes()):
        G.nodes[i]['demand']    = pDemand['MW'][i]
        G.nodes[i]['ColorNode'] = '#5890BD'
    
    # nx.layout.shell_layout only works for more than 3 nodes
    if len(shell2)>1:
#         pos = nx.drawing.layout.shell_layout(G, shells)
        pos = {i: (dfNodeLocation['Latitude'][i], dfNodeLocation['Longitude'][i]) for i in list(G.nodes())}
    else:
        pos = nx.drawing.layout.spring_layout(G)
    for node in G.nodes:
        G.nodes[node]['pos'] = list(pos[node])


    if len(shell2)==0:
        traceRecode = []  # contains edge_trace, node_trace, middle_node_trace

        node_trace = go.Scatter(x=tuple([1]), y=tuple([1]), text=tuple([str(InputNode)]), textposition="bottom center",
                                mode='markers+text',
                                marker={'size': 10, 'color': 'LightSkyBlue'})
        traceRecode.append(node_trace)

        node_trace1 = go.Scatter(x=tuple([1]), y=tuple([1]),
                                mode='markers',
                                marker={'size': 10, 'color': 'LightSkyBlue'},
                                opacity=0)
        traceRecode.append(node_trace1)

        figure = {
            "data": traceRecode,
            "layout": go.Layout(title='Interactive Visualization', showlegend=False,
                                margin={'b': 40, 'l': 40, 'r': 40, 't': 40},
                                xaxis={'showgrid': False, 'zeroline': False, 'showticklabels': False},
                                yaxis={'showgrid': False, 'zeroline': False, 'showticklabels': False},
                                height=600
                                )}
        return figure
    traceRecode = []  # contains edge_trace, node_trace, middle_node_trace
############################################################################################################################################################
    index = 0
    for edge in G.edges:
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
#         weight = float(G.edges[edge]['MW']) / max(dfNetwork['MW']) * 10
        weight = float(G.edges[edge]['MW']) * 10
        trace = go.Scatter(x=tuple([x0, x1, None]), y=tuple([y0, y1, None]),
                           mode='lines',
                           line={'width': weight},
                           marker=dict(color=G.edges[edge]['ColorType']),
                           line_shape='spline',
                           opacity=1)
        traceRecode.append(trace)
        index = index + 1
###############################################################################################################################################################
    node_trace = go.Scatter(x=[], y=[], hovertext=[], text=[], mode='markers+text', textposition="bottom center",
                            hoverinfo="text", marker={'size': 10, 'color': 'LightSkyBlue'})
    index = 0
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        hovertext = "GenOut: " + str(round(dfGenToNode['MW'][node],2)) + " MW" + "<br>" + "Demand: " + str(round(pDemand['MW'][node],2)) + " MW" + "<br>" + "PNS: " + str(round(pPNS['MW'][node],2)) + " MW"
        text = node
        node_trace['x'] += tuple([x])
        node_trace['y'] += tuple([y])
        node_trace['hovertext'] += tuple([hovertext])
        node_trace['text'] += tuple([text])
        index = index + 1
    traceRecode.append(node_trace)
################################################################################################################################################################
    middle_hover_trace = go.Scatter(x=[], y=[], hovertext=[], mode='markers', hoverinfo="text",
                                    marker={'size': 10, 'color': 'LightSkyBlue'},
                                    opacity=0)
    index = 0
    for edge in G.edges:
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
        hovertext = "From: " + str(G.edges[edge]['InitialNode']) + "<br>" + "To: " + str(G.edges[edge]['FinalNode']) + "<br>" + "NetworkFlow: " + str(round(G.edges[edge]['vFlow'],2))+ " MW" + "<br>" + "NetworkUse: " + str(round(G.edges[edge]['MW']*100,2)) + " %"
        middle_hover_trace['x'] += tuple([(x0 + x1) / 2])
        middle_hover_trace['y'] += tuple([(y0 + y1) / 2])
        middle_hover_trace['hovertext'] += tuple([hovertext])
        index = index + 1
    traceRecode.append(middle_hover_trace)
#################################################################################################################################################################
    figure = {
        "data": traceRecode,
        "layout": go.Layout(title='Interactive Visualization -- Date: '+ str(InputDate), showlegend=False, hovermode='closest',
                            margin={'b': 40, 'l': 40, 'r': 40, 't': 40},
                            xaxis={'showgrid': False, 'zeroline': False, 'showticklabels': False},
                            yaxis={'showgrid': False, 'zeroline': False, 'showticklabels': False},
                            height=600,
                            clickmode='event+select',
                            annotations=[
                                dict(
                                    ax=(G.nodes[edge[0]]['pos'][0] + G.nodes[edge[1]]['pos'][0]) / 2,
                                    ay=(G.nodes[edge[0]]['pos'][1] + G.nodes[edge[1]]['pos'][1]) / 2, axref='x', ayref='y',
                                    x=(G.nodes[edge[1]]['pos'][0] * 3 + G.nodes[edge[0]]['pos'][0]) / 4,
                                    y=(G.nodes[edge[1]]['pos'][1] * 3 + G.nodes[edge[0]]['pos'][1]) / 4, xref='x', yref='y',
                                    showarrow=True,
                                    arrowhead=3,
                                    arrowsize=4,
                                    arrowwidth=1,
                                    opacity=1
                                ) for edge in G.edges]
                            )}
    return figure
######################################################################################################################################################################
#Date Range Selector
dfNetworkDate              = pd.read_csv('oT_Result_NetworkUtilization_ToPT_' + CaseName + '.csv')
dfNetworkDate['LoadLevel'] = pd.to_datetime(dfNetworkDate.LoadLevel)
######################################################################################################################################################################
# styles: for right side hover/click component
styles = {
    'pre': {
        'border': 'thin lightgrey solid',
        'overflowX': 'scroll'
    }
}

app.layout = html.Div([
    #########################Title
    html.Div([html.H1(CaseName+": "+"Network Graph")],
             className="row",
             style={'textAlign': "center"}),
    #############################################################################################define the row
    html.Div(
        className="row",
        children=[
            ##############################################left side two input components
            html.Div(
                className="two columns",
                children=[
                    html.Div([
                        dcc.Markdown(d("""
                        **Select the date to visualize...**
                        """)),
                        dcc.DatePickerSingle(
                          id='my-date-picker-single',
                          min_date_allowed     =datetime(dfNetworkDate['LoadLevel'].max().to_pydatetime().year, 1, 1),
                          max_date_allowed     =datetime(dfNetworkDate['LoadLevel'].max().to_pydatetime().year,
                                                         dfNetworkDate['LoadLevel'].max().to_pydatetime().month,
                                                         dfNetworkDate['LoadLevel'].max().to_pydatetime().day),
                          initial_visible_month=datetime(dfNetworkDate['LoadLevel'].max().to_pydatetime().year,
                                                         dfNetworkDate['LoadLevel'].max().to_pydatetime().month, 1),
                          date                 =str(datetime(dfNetworkDate['LoadLevel'].max().to_pydatetime().year,
                                                         dfNetworkDate['LoadLevel'].max().to_pydatetime().month,
                                                         dfNetworkDate['LoadLevel'].max().to_pydatetime().day, 23, 59, 59)),
                        ),
                        html.Div(id='output-container-date-picker-single', style={'display': 'none'})
                        ], className="row ", style={'marginTop': 30, 'marginBottom': 15}),
                    html.Div([html.Button('Submit', id='btn-nclicks-1', n_clicks=0),
                             html.Div(id='container-button-timestamp', style={'display': 'none'})])
                ]
            ),

            ############################################middle graph component
            html.Div(
                className="eight columns",
                children=[
                    dcc.Markdown(d("""
                            **Slide the bar to define the hour...**
                            """)),
                    # Radio Button
                    dcc.Slider(
                        id='my-slider',
                        min=0,
                        max=23,
                        step=None,
                        marks={
                            0: '00:00',
                            1: '01:00',
                            2: '02:00',
                            3: '03:00',
                            4: '04:00',
                            5: '05:00',
                            6: '06:00',
                            7: '07:00',
                            8: '08:00',
                            9: '09:00',
                            10: '10:00',
                            11: '11:00',
                            12: '12:00',
                            13: '13:00',
                            14: '14:00',
                            15: '15:00',
                            16: '16:00',
                            17: '17:00',
                            18: '18:00',
                            19: '19:00',
                            20: '20:00',
                            21: '21:00',
                            22: '22:00',
                            23: '23:00'
                        },
                        value=5
                    ),
                    html.Div(id='slider-output-container', style={'display': 'none'}),
                    dcc.Graph(id="my-graph",
                                    figure=network_graph(InputDate, InputTime))],
            ),

            #########################################right side two output component
            html.Div(
                className="two columns",
                children=[
                    html.Div(
                        className='twelve columns',
                        children=[
                            dcc.Markdown(d("""
                            **Hover Data**
                            Mouse over values in the graph.
                            """)),
                            html.Pre(id='hover-data', style=styles['pre'])
                        ],
                        style={'height': '400px'}),

                    html.Div(
                        className='twelve columns',
                        children=[
                            dcc.Markdown(d("""
                            **Click Data**
                            Click on points in the graph.
                            """)),
                            html.Pre(id='click-data', style=styles['pre'])
                        ],
                        style={'height': '400px'})
                ]
            )
        ]
    )
])

###################################callback for left side components


@app.callback(
    dash.dependencies.Output('slider-output-container', 'children'),
    [dash.dependencies.Input('my-slider', 'value')])
def update_output(value):
    return value

@app.callback(
    dash.dependencies.Output('output-container-date-picker-single', 'children'),
    [dash.dependencies.Input('my-date-picker-single', 'date')])
def update_output(date):
    if date is not None:
        date = datetime.strptime(re.split('T| ', date)[0], '%Y-%m-%d')
        date = datetime(date.year, date.month, date.day, 0, 0, 0)
        return date

@app.callback(dash.dependencies.Output('container-button-timestamp', 'children'),
              [dash.dependencies.Input('btn-nclicks-1', 'n_clicks'),
               dash.dependencies.Input('slider-output-container', 'children'),
               dash.dependencies.Input('output-container-date-picker-single', 'children')])
def displayClick(btn1, InputValue, InputDate):
    changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
    if 'btn-nclicks-1' in changed_id:
        msgValue = InputValue
        msgDate  = datetime.strptime(re.split('T| ', InputDate)[0], '%Y-%m-%d')
        msgDate  = datetime(msgDate.year, msgDate.month, msgDate.day, 0, 0, 0)
        msg      = [msgDate, msgValue]
    else:
        msg = ["2030-03-31 01:00:00", 0]
    return msg

@app.callback(
    dash.dependencies.Output('my-graph', 'figure'),
    [dash.dependencies.Input('container-button-timestamp', 'children')])
def update_output(Inputmsg):
    return network_graph(Inputmsg[0], Inputmsg[1])
#     to update the global variable of YEAR and ACCOUNT
    
################################callback for right side components
@app.callback(
    dash.dependencies.Output('hover-data', 'children'),
    [dash.dependencies.Input('my-graph', 'hoverData')])
def display_hover_data(hoverData):
    return json.dumps(hoverData, indent=2)


@app.callback(
    dash.dependencies.Output('click-data', 'children'),
    [dash.dependencies.Input('my-graph', 'clickData')])
def display_click_data(clickData):
    return json.dumps(clickData, indent=2)



if __name__ == '__main__':
    app.run_server(debug=True)
    
# if __name__ == '__main__':
# #     app.run_server(debug=True)
#     app.server.run(port=8000, host='127.0.0.1')