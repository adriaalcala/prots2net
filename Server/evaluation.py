import csv
import math


def compute_confusion_matrix(predicted_edges, real_edges, confidence_threshold, score_threshold, num_prots):
    predicted = {tuple(sorted((e['protein_1'], e['protein_2']))) for e in predicted_edges if float(e['confidence']) >= confidence_threshold}
    real = {tuple(sorted((e['protein1'], e['protein2']))) for e in real_edges if float(e['combined_score']) >= score_threshold}
    TP = len(predicted.intersection(real))
    FP = len(predicted.difference(real))
    FN = len(real.difference(predicted))
    TN = ((num_prots * (num_prots - 1)) / 2) - FN
    return TP, FP, FN, TN

def roc_curve(real_file, predicted_file, score_threshold=700):
    with open(real_file) as f:
        reader = csv.DictReader(f, delimiter=' ')
        real = list(reader)
    with open(predicted_file) as f:
        reader = csv.DictReader(f, delimiter=' ')
        predicted = list(reader)
    num_prots = len(set(e['protein1'] for e in real).union(set(e['protein2'] for e in real)))
    return [compute_confusion_matrix(predicted, real, confidence / 100, score_threshold, num_prots) for confidence in range(101)]


def accuracy(roc_curve):
    return [(r[0] + r[3])/(sum(r)) for r in roc_curve]

def sensitivity(roc_curve):
    return [r[0]/(r[0] + r[2]) for r in roc_curve]

def precision(roc_curve):
    return [r[0]/(r[0] + r[1]) if r[0] else 0 for r in roc_curve]

def mcc(roc_curve):
    return [((r[0] * r[3]) - (r[1] * r[2]))/(math.sqrt((r[0] + r[1]) * (r[0] + r[2]) * (r[3] + r[1]) * (r[3] + r[2]))) if (r[0] + r[1]) * (r[0] + r[2]) * (r[3] + r[1]) * (r[3] + r[2]) else 0 for r in roc_curve]