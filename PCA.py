import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.metrics import silhouette_score
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

from tqdm.auto import tqdm
from PIL import Image
import os

data = pd.read_csv('data/orf.csv')
data = data.dropna()

# Build distance matrix
array = []
listofpoints = []
for experiment in data['step'].unique():
    if experiment == 'P30': continue
    for replicate in data['replicate'].unique():
        orfs = np.zeros(156)
        for orf in data[(data['step'] == experiment) & (data['replicate'] == replicate)]['orf'].values:
            orfs[int(orf.replace("ORF", ""))-1] = data[(data['step'] == experiment) &
                                                       (data['replicate'] == replicate) &
                                                       (data['orf'] == orf)]['count'].values[0]
        orfs = orfs / orfs.sum()
        orfs = np.nan_to_num(orfs)
        array.append(orfs)
        listofpoints.append((experiment, int(replicate)))

# {P15, P30, P50, P65, P95} = 5
# {1, 2, 3, 4, 5, 6, 7, 8, 9, 10} = 10
# 5 * 10 = 50
# On fera le tri chaud froid après

# Distance matrix
size = len(array)
diss_matrix = np.zeros((size, size))
for i in range(size):
    for j in range(size):
        diss_matrix[i, j] = np.linalg.norm(array[i] - array[j])

# PCA
pca = PCA(n_components=3)
pca_result = pca.fit_transform(array)

# Plot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

forms = ['o', 's', 'D', 'x', '+']
temp = ['blue', 'red']
steps = ['P15', 'P50', 'P65', 'P90']

for i, (x, y, z) in enumerate(pca_result):
    marker = forms[i // 10]
    c = temp[i % 2] if i > 9 else 'green'
    label = f"{steps[i // 10]}" if i // 10 != (i-1) // 10 else None
    ax.scatter(x, y, z, c=c, marker=marker, label=label)

ax.set_xlabel(f'PC1 {pca.explained_variance_ratio_[0] * 100:.2f}%')
ax.set_ylabel(f'PC2 {pca.explained_variance_ratio_[1] * 100:.2f}%')
ax.set_zlabel(f'PC3 {pca.explained_variance_ratio_[2] * 100:.2f}%')
# ax.set_title(f'PCA {pca.explained_variance_ratio_.sum() * 100:.2f}%')
print(f'PCA {pca.explained_variance_ratio_.sum() * 100:.2f}%')
plt.legend()

def generate_frames(path='test'):
    if not os.path.exists(path): os.makedirs(path)
    for angle in tqdm(range(0, 360*3 + 1), desc='Generating frames'):
        # Normalize the angle to the range [-180, 180] for display
        angle_norm = (angle + 180) % 360 - 180

        # Cycle through a full rotation of elevation, then azimuth, roll, and all
        elev = azim = roll = 0
        if angle <= 360:
            elev = angle_norm
        elif angle <= 360*2:
            azim = angle_norm
        else:
            elev = azim = roll = angle_norm

        # Update the axis view and title
        ax.view_init(elev, azim, roll)
        # plt.draw()
        plt.pause(0.01)
        fig.tight_layout()
        if angle < 10:
            fig.savefig(f'{path}/000{angle}.png', dpi=300)
        elif angle < 100:
            fig.savefig(f'{path}/00{angle}.png', dpi=300)
        elif angle < 1000:
            fig.savefig(f'{path}/0{angle}.png', dpi=300)
        else:
            fig.savefig(f'{path}/{angle}.png', dpi=300)

def make_gif(dirpath='PCA'):
    frames = []
    paths = os.listdir(dirpath)
    paths.sort()
    # take half the frames
    paths = paths[::2]
    for path in tqdm(paths, desc='Loading frames'):
        img = Image.open(dirpath + "/" + path)
        img = img.resize((640, 480))
        frames.append(img)

    starter = frames.pop(0)
    # make it loop
    starter.save(f"{dirpath}.gif", save_all=True, append_images=frames, loop=0)
    # save as mp4
    os.system(f'ffmpeg -y -i {dirpath}.gif -vf "scale=640:-1:flags=lanczos" {dirpath}.mp4')

### Uncomment to generate the frames and make the gif ###
#generate_frames('PCA')
#make_gif('PCA')
#########################################################

plt.close()
def gap_statistic(X, max_k=15, n_references=10):
    gaps = []
    for k in range(1, max_k):
        # Fit KMeans to real data
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10).fit(X)
        Wk = np.log(sum(np.min(cdist(X, kmeans.cluster_centers_, 'euclidean'), axis=1)) / X.shape[0])

        # Generate reference datasets
        Wk_refs = []
        for _ in range(n_references):
            random_data = np.random.uniform(X.min(axis=0), X.max(axis=0), X.shape)
            kmeans_ref = KMeans(n_clusters=k, random_state=42, n_init=10).fit(random_data)
            Wk_refs.append(
                np.log(sum(np.min(cdist(random_data, kmeans_ref.cluster_centers_, 'euclidean'), axis=1)) / X.shape[0]))

        gap = np.mean(Wk_refs) - Wk
        gaps.append(gap)

    return gaps


# Compute the gap statistic for PCA-transformed data
gaps = gap_statistic(pca_result)

# Plot the Gap Statistic
plt.plot(range(1, len(gaps) + 1), gaps, marker='o')
plt.xlabel('Number of Clusters (K)')
plt.ylabel('Gap Statistic')
plt.title('Gap Statistic for Optimal K')
plt.tight_layout()
plt.show()
plt.close()

silhouette = []
for k in range(2, 15):  # Focus on reasonable K values
    kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
    labels = kmeans.fit_predict(pca_result)
    score = silhouette_score(pca_result, labels)
    silhouette.append(score)

plt.plot(range(2, 15), silhouette, marker='o')
plt.xlabel('Number of Clusters (K)')
plt.ylabel('Silhouette Score')
plt.title('Silhouette Score for Optimal K')
plt.tight_layout()
plt.show()
plt.close()

n_clusters = 2

# Apply K-Means
kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
clusters = kmeans.fit_predict(pca_result)

p_to_c = {n: [] for n in range(n_clusters)}
for i, c in enumerate(clusters):
    p_to_c[int(c)].append(listofpoints[i])


# 3D Scatter Plot
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Assign colors to clusters
colors = ['orange', 'cyan', 'purple', 'green', 'blue', 'red', 'black', 'yellow', 'pink', 'brown']
for i in range(n_clusters):  # Loop through clusters
    ax.scatter(
        pca_result[clusters == i, 0],
        pca_result[clusters == i, 1],
        pca_result[clusters == i, 2],
        c=colors[i],
        label=f'Cluster {i+1}',
        alpha=0.7
    )

# Plot cluster centers
ax.scatter(kmeans.cluster_centers_[:, 0],
           kmeans.cluster_centers_[:, 1],
           kmeans.cluster_centers_[:, 2],
           c='black', marker='X', s=200, label='Barycentre')

# Labels
ax.set_xlabel(f'PC1 {pca.explained_variance_ratio_[0] * 100:.2f}%')
ax.set_ylabel(f'PC2 {pca.explained_variance_ratio_[1] * 100:.2f}%')
ax.set_zlabel(f'PC3 {pca.explained_variance_ratio_[2] * 100:.2f}%')

plt.tight_layout()
plt.legend()

### Uncomment to generate the frames and make the gif ###
#generate_frames('kmeans')
#make_gif('kmeans')
#########################################################
