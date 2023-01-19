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
