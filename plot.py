import matplotlib.pyplot as plt
from collections import Counter

def plot_codon_repartition(bact_name, start, stop):

    c_start = Counter(start)
    c_stop = Counter(stop)

    #Piechart start
    plt.pie(c_start.values(), labels=c_start.keys(), labeldistance=1.15)
    plt.title(f"Start codon repartition of {bact_name}")
    plt.show()
    #Piechart stop
    plt.pie(c_stop.values(), labels=c_stop.keys(), labeldistance=1.15)#modif les label : ajouter catégorie autres
    plt.title(f"Stop codon repartition of {bact_name}")
    plt.semilogy() #ajouter une échelle log pour la visibilité
    plt.show()

    #Barplot start
    y_pos = range(len(c_start.keys()))
    plt.bar(y_pos, c_start.values())
    plt.xticks(y_pos, c_start.keys(),fontsize=5,rotation=50)
    plt.title(f"Start codon repartition of {bact_name}")
    plt.show()
    #Barplot stop
    y_pos = range(len(c_stop.keys()))
    plt.bar(y_pos, c_stop.values())
    plt.xticks(y_pos, c_stop.keys(),fontsize=5,rotation=50)
    plt.title(f"Stop codon repartition of {bact_name}")
    plt.semilogy()
    plt.show()  