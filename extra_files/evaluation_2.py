import csv
import math
import random
from sklearn import metrics
import plotly.express as px
import plotly.figure_factory as ff

import pandas as pd


def compute_confusion_matrix(predicted_edges, real, negative, confidence_threshold):
    if confidence_threshold == 0:
        return len(real), len(negative), 0, 0
    predicted = {tuple(sorted((e['protein_1'], e['protein_2']))) for e in predicted_edges if float(e['confidence']) >= confidence_threshold}
    TP = len(predicted.intersection(real))
    FP = len(predicted.intersection(negative))
    FN = len(real.difference(predicted))
    TN = len(negative) - FP
    return TP, FP, FN, TN


def roc_curve(real_file, predicted_file, score_threshold=999, not_sure_threshold=200, percent_positives=1):
    with open(real_file) as f:
        reader = csv.DictReader(f, delimiter=' ')
        real = list(reader)
    with open(predicted_file) as f:
        reader = csv.DictReader(f)
        predicted = list(reader)
    prots = set(e['protein1'] for e in real).union(set(e['protein2'] for e in real))
    sure_real = {tuple(sorted((e['protein1'], e['protein2']))) for e in real if float(e['combined_score']) >= score_threshold}
    not_sure = {tuple(sorted((e['protein1'], e['protein2']))) for e in real if float(e['combined_score']) >= not_sure_threshold}
    # real = set(random.sample(real, int(len(real) * percent_positives)))
    # negatives = {tuple(sorted((p1, p2))) for p1 in prots for p2 in prots if p1 != p2 and tuple(sorted((p1, p2))) not in real}
    # negatives = set(random.sample(negatives, len(real)))
    negatives = set()
    num_edges = len(sure_real)
    while len(negatives) < num_edges:
        p1, p2 = random.sample(prots, 2)
        if tuple(sorted((p1, p2))) not in not_sure:
            negatives.add(tuple(sorted((p1, p2))))
        # print(len(negatives), num_edges)
    return [compute_confusion_matrix(predicted, sure_real, negatives, confidence / 1000,) for confidence in range(1001)]

def metrics(real_file, predicted_file, score_threshold=999, not_sure_threshold=200, percent_positives=1):
    with open(real_file) as f:
        reader = csv.DictReader(f, delimiter=' ')
        real = list(reader)
    with open(predicted_file) as f:
        reader = csv.DictReader(f)
        predicted = list(reader)
    prots = set(e['protein1'] for e in real).union(set(e['protein2'] for e in real))
    sure_real = {tuple(sorted((e['protein1'], e['protein2']))) for e in real if float(e['combined_score']) >= score_threshold}
    not_sure = {tuple(sorted((e['protein1'], e['protein2']))) for e in real if float(e['combined_score']) >= not_sure_threshold}
    # real = set(random.sample(real, int(len(real) * percent_positives)))
    # negatives = {tuple(sorted((p1, p2))) for p1 in prots for p2 in prots if p1 != p2 and tuple(sorted((p1, p2))) not in real}
    # negatives = set(random.sample(negatives, len(real)))
    negatives = set()
    num_edges = len(sure_real)
    while len(negatives) < num_edges:
        p1, p2 = random.sample(prots, 2)
        if tuple(sorted((p1, p2))) not in not_sure:
            negatives.add(tuple(sorted((p1, p2))))
        # print(len(negatives), num_edges)
    r = [compute_confusion_matrix(predicted, sure_real, negatives, 0.5),]
    return accuracy(r), sensitivity(r), specificity(r), precision(r), mcc(r), fpr(r), f1(r)


def accuracy(roc_curve):
    return [(r[0] + r[3])/(sum(r)) for r in roc_curve]

def sensitivity(roc_curve):
    return [r[0]/(r[0] + r[2]) for r in roc_curve]

def specificity(roc_curve):
    return [r[3]/(r[1] + r[3]) for r in roc_curve]


def precision(roc_curve):
    return [r[0]/(r[0] + r[1]) if r[0] else 0 for r in roc_curve]


def mcc(roc_curve):
    return [((r[0] * r[3]) - (r[1] * r[2]))/(math.sqrt((r[0] + r[1]) * (r[0] + r[2]) * (r[3] + r[1]) * (r[3] + r[2]))) if (r[0] + r[1]) * (r[0] + r[2]) * (r[3] + r[1]) * (r[3] + r[2]) else 0 for r in roc_curve]

def fpr(roc_curve):
    return [r[1]/(r[1] + r[3]) for r in roc_curve]

def f1(roc_curve):
    return [(2 * r[0])/(2*r[0] + r[1] + r[2]) for r in roc_curve]


def plot_confusion_matrix(r):
    print(r)
    z = [[r[3], r[1]], [r[2], r[0]]]
    t = r[0] + r[2]
    print(r, t, z)
    z_2 = [[100 * z[0][0]/t, 100 * z[0][1]/t], [100 * z[1][0]/t, 100 * z[1][1]/t]]
    x = y = ['No Interact', 'Interact']
    z_text = [[f"{z[0][0]} ({z_2[0][0]:.2f}%)", f"{z[0][1]} ({z_2[0][1]:.2f}%)"], [f"{z[1][0]} ({z_2[1][0]:.2f}%)", f"{z[1][1]} ({z_2[1][1]:.2f}%)"]]
    fig = ff.create_annotated_heatmap(z, x=x, y=y, annotation_text=z_text, colorscale='BuGn')
    # add title
    fig.update_layout(title_text='<i><b>Confusion matrix</b></i>') #, xaxis = dict(title='Predicted value'), yaxis = dict(title='Real value'))
    # add custom xaxis title
    fig.add_annotation(dict(font=dict(color="black",size=14), x=0.5, y=-0.15, showarrow=False, text="Predicted value", xref="paper", yref="paper"))
    # add custom yaxis title
    fig.add_annotation(dict(font=dict(color="black",size=14), x=-0.35, y=0.5, showarrow=False, text="Real value", textangle=-90, xref="paper", yref="paper"))
    # adjust margins to make room for yaxis title
    fig.update_layout(margin=dict(t=50, l=200))
    # add colorbar
    fig['data'][0]['showscale'] = True
    fig.write_image(f"images/confusion_matrix_{species_1}_{species_2}.png")

def create_figures(species_1, species_2):
    r = roc_curve(f'data/edges/{species_1}.protein.links.v11.5.txt', f'data/nets/Net_{species_1}_{species_2}.csv')
    threshold = [i/1000 for i in range(1001)]
    a = accuracy(r)
    df_a = pd.DataFrame({'threshold': threshold[1:-1], 'score': a[1:-1], 'metric': 'accuracy'})
    p = precision(r)
    df_p = pd.DataFrame({'threshold': threshold[1:-1], 'score': p[1:-1], 'metric': 'precision'})
    s = sensitivity(r)
    df_s = pd.DataFrame({'threshold': threshold[1:-1], 'score': s[1:-1], 'metric': 'sensitivity'})
    sp = specificity(r)
    df_sp = pd.DataFrame({'threshold': threshold[1:-1], 'score': sp[1:-1], 'metric': 'specificity'})
    mcc_values = mcc(r)
    df_mcc = pd.DataFrame({'threshold': threshold[1:-1], 'score': mcc_values[1:-1], 'metric': 'MCC'})
    f1_values = f1(r)
    df_f1 = pd.DataFrame({'threshold': threshold[1:-1], 'score': f1_values[1:-1], 'metric': 'F1-score'})
    f = fpr(r)
    df_f1 = pd.DataFrame({'threshold': threshold[1:-1], 'score': f[1:-1], 'metric': 'FallOut'})

    df = pd.concat([df_a, df_p, df_s, df_sp, df_mcc, df_f1])

    fig = px.line(df, x='threshold', y='score', color='metric')
    fig.write_image(f"images/metrics_{species_1}_{species_2}.png")
    t = sensitivity(r)
    f2 = [i + (101-idx) * 1e-5 for idx, i in enumerate(f)]
    auc = metrics.auc(f2, t)
    print(auc)
    fig = px.area(
        x=f2, y=t,
        title=f'ROC Curve (AUC={auc:.4f})',
        labels=dict(x='False Positive Rate', y='True Positive Rate'),
        width=700, height=500
    )
    fig.add_shape(
        type='line', line=dict(dash='dash'),
        x0=0, x1=1, y0=0, y1=1
    )

    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    fig.update_xaxes(constrain='domain')
    fig.write_image(f"images/ROC_{species_1}_{species_2}.png")
    s2 = [i + (101-idx) * 1e-5 for idx, i in enumerate(s)]
    p2 = [i - idx * 1e-5 for idx, i in enumerate(p)]
    auc_2 = metrics.auc(s2, p2)
    fig = px.area(
        x=s2[1:-1], y=p[1:-1],
        title=f'Precision-Recall Curve (AUC={auc_2:.4f})',
        labels=dict(x='Recall', y='Precision'),
        width=700, height=500
    )
    fig.add_shape(
        type='line', line=dict(dash='dash'),
        x0=0, x1=1, y0=0, y1=1
    )

    fig.update_yaxes(scaleanchor="x", scaleratio=1)
    fig.update_xaxes(constrain='domain')
    fig.write_image(f"images/precision_recall_{species_1}_{species_2}.png")

    plot_confusion_matrix(r[500])

if __name__ == '__main__':
    species = [(9606, 3702), (9606, 224325), [9606, 9597], [4932, 4920], [633813, 309807], [9606, 4920], [9606, 309807]]
    for species_1, species_2 in species:
        print(species_1, species_2)
        print(metrics(f'data/edges/{species_1}.protein.links.v11.5.txt', f'data/nets/Net_{species_1}_{species_2}.csv'))

