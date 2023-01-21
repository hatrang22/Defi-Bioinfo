from sklearn.cluster import KMeans, AgglomerativeClustering
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import RidgeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score

def clustering(dfs, n_cluster, method, top3_only=False):
    
    top3 = ["ATG", "GTG", "TTG"]

    method = method.lower()

    df = pd.concat(dfs, ignore_index=True).drop("ID", axis=1)
    
    if top3_only:
        df = df[top3]

    df = df.fillna(0)
    
    X = df.to_numpy()

    if method == "kmeans":
        opt = KMeans(n_clusters=n_cluster)

    elif method == "agglomerativeclustering":
        opt = AgglomerativeClustering(n_clusters=n_cluster, compute_distances=True)

    else:
        raise ValueError(f"Unknown clustering method {method}")

    opt = opt.fit(X)

    return opt

def classify(dfs, method, test_size=0.33, split_seed=0, training_seed=0, top3_only=False):

    top3 = ["ATG", "GTG", "TTG"]

    method = method.lower()

    df = pd.concat(dfs, ignore_index=True).drop("ID", axis=1)
    
    if top3_only:
        df = df[top3]

    df = df.fillna(0)
    
    X = df.to_numpy()

    y = [] # create target label
    for i in range(len(dfs)):
        y += [i for j in range(len(dfs[i]))]
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=split_seed)
    
    if method == "linear":
        model = RidgeClassifier(random_state=training_seed)
    
    elif method == "decisiontree":
        model = DecisionTreeClassifier(random_state=training_seed)
    
    elif method == "randomforest":
        model = RandomForestClassifier(random_state=training_seed)

    elif method == "kneighbors":
        model = KNeighborsClassifier()

    elif method == "svm":
        model = SVC(random_state=training_seed)

    elif method == "neuralnetwork":
        model = MLPClassifier(hidden_layer_sizes=(X.shape[-1],),  # number of neurons based on number of features 
                                solver="lbfgs",  # for small datasets, lbfgs solver can work better
                                early_stopping=True,  # to avoid overfitting using validation set
                                random_state=training_seed,  # random seed to initialize weights
                                max_iter=400)  # max iteration
    
    else:
        raise ValueError(f"Unknown classification method {method}")

    model = model.fit(X_train, y_train)

    y_pred = model.predict(X_test)

    mat = confusion_matrix(y_test, y_pred)

    global_score = accuracy_score(y_test, y_pred)

    return mat, global_score

#%% ROC curve
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt

def calculate_tpr_fpr(y_real, y_pred):
    '''
    Calculates the True Positive Rate (tpr) and the True Negative Rate (fpr) based on real and predicted observations
    
    Args:
        y_real: The list or series with the real classes
        y_pred: The list or series with the predicted classes
        
    Returns:
        tpr: The True Positive Rate of the classifier
        fpr: The False Positive Rate of the classifier
    '''
    
    # Calculates the confusion matrix and recover each element
    cm = confusion_matrix(y_real, y_pred)
    TN = cm[0, 0]
    FP = cm[0, 1]
    FN = cm[1, 0]
    TP = cm[1, 1]
    
    # Calculates tpr and fpr
    tpr =  TP/(TP + FN) # sensitivity - true positive rate
    fpr = 1 - TN/(TN+FP) # 1-specificity - false positive rate
    
    return tpr, fpr

def get_all_roc_coordinates(y_real, y_proba):
    '''
    Calculates all the ROC Curve coordinates (tpr and fpr) by considering each point as a threshold for the predicion of the class.
    
    Args:
        y_real: The list or series with the real classes.
        y_proba: The array with the probabilities for each class, obtained by using the `.predict_proba()` method.
        
    Returns:
        tpr_list: The list of TPRs representing each threshold.
        fpr_list: The list of FPRs representing each threshold.
    '''
    tpr_list = [0]
    fpr_list = [0]
    for i in range(len(y_proba)):
        threshold = y_proba[i]
        y_pred = y_proba >= threshold
        tpr, fpr = calculate_tpr_fpr(y_real, y_pred)
        tpr_list.append(tpr)
        fpr_list.append(fpr)
    return tpr_list, fpr_list


def plot_roc_curve(tpr, fpr, classes, scatter = False):
    '''
    Plots the ROC Curve by using the list of coordinates (tpr and fpr).
    
    Args:
        tpr: The list of TPRs representing each coordinate.
        fpr: The list of FPRs representing each coordinate.
        scatter: When True, the points used on the calculation will be plotted with the line (default = True).
    '''
    if scatter:
        plt.scatter(x = fpr, y = tpr)
    
    fpr_tpr = pd.DataFrame({'fpr': fpr, 'tpr': tpr})
    fpr_tpr.sort_values('fpr', inplace=True)

    plt.plot(fpr_tpr['fpr'], fpr_tpr['tpr'])
    plt.plot([0,1],[0,1],'r') #diagonal line
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title('ROC curve of '+ classes)
    plt.show()
    
def ROC_curves_and_AOC_scores(dfs):
    fusion=pd.concat(dfs, ignore_index=True)
    X_=fusion[['ATC','GTG','TTG','CTG']]
    y_=fusion['phylum']
    #split jeu de données
    x_train, x_test, y_train, y_test  = train_test_split(X_, 
                                                         y_, 
                                                         test_size=0.25, 
                                                         random_state=42)
    #création du modèle
    modele_rf = KNeighborsClassifier()
    #apprentissage
    modele_rf.fit(x_train, y_train)
    #predict proba
    y_scores = modele_rf.predict_proba(x_test)
    
    #ROC plot and ROC_auc_scores print
    classes = modele_rf.classes_
    roc_auc_ovr = {}

    for i in range(len(classes)):
        # Gets the class
        c = classes[i]
        
        # Prepares an auxiliar dataframe to help with the plots
        df_aux = x_test.copy()
        df_aux['class'] = [1 if y == c else 0 for y in y_test]
        df_aux['prob'] = y_scores[:, i]
        df_aux = df_aux.reset_index(drop = True)
        
        # Calculates the ROC Coordinates and plots the ROC Curves
        tpr, fpr = get_all_roc_coordinates(df_aux['class'], df_aux['prob'])
        plot_roc_curve(tpr, fpr, scatter = False, classes=c)
        
        # Calculates the ROC AUC OvR
        roc_auc_ovr[c] = roc_auc_score(df_aux['class'], df_aux['prob'])

    # Displays the ROC AUC for each class
    avg_roc_auc = 0
    i = 0
    for k in roc_auc_ovr:
        avg_roc_auc += roc_auc_ovr[k]
        i += 1
        print(f"{k} ROC AUC OvR: {roc_auc_ovr[k]:.4f}")
    print(f"average ROC AUC OvR: {avg_roc_auc/i:.4f}")
