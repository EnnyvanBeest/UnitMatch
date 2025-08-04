from sklearn.neighbors import KernelDensity
import numpy as np
from tqdm import tqdm


class NaiveBayes:
    def __init__(self, distance_bandwidth=50):
        self.distance_bandwidth = distance_bandwidth
        # Similarity will use default bandwidth from sklearn
        self.class_priors = None
        self.similarity_kde = None
        self.distance_kde = None
        
    def fit(self, X, y):
        # X[:, 0] should be similarity
        # X[:, 1] should be distance
        
        # Calculate class priors
        classes = np.unique(y)
        n_samples = len(y)
        self.classes = classes
        self.class_priors = {}

        same_units = X[y == True]
        self.class_priors[True] = len(same_units) / n_samples

        similarity_values = same_units[:, 0].reshape(-1, 1)
        distance_values = same_units[:, 1].reshape(-1, 1)

        # Fit KDE for similarity (default bandwidth)
        sim_kde = KernelDensity(kernel='gaussian')
        sim_kde.fit(similarity_values)
        self.similarity_kde = sim_kde
        
        # Fit KDE for distance (bandwidth specified in constructor)
        dist_kde = KernelDensity(kernel='gaussian', bandwidth=self.distance_bandwidth)
        dist_kde.fit(distance_values)
        self.distance_kde = dist_kde
    
    def predict_prob(self, X):
        n_samples = X.shape[0]
        probs = np.zeros((n_samples,))
        
        # For each sample
        for i, sample in tqdm(enumerate(X), total=X.shape[0]):
            similarity = sample[0].reshape(-1, 1)
            distance = sample[1].reshape(-1, 1)

            # Start with log of prior probability
            log_likelihood = np.log(self.class_priors[True])
            
            # Add log probability from similarity KDE
            sim_log_prob = self.similarity_kde.score_samples(similarity)[0]
            log_likelihood += sim_log_prob
            
            # Add log probability from distance KDE
            dist_log_prob = self.distance_kde.score_samples(distance)[0]
            log_likelihood += dist_log_prob
            
            probs[i] = log_likelihood  # Still in log space to avoid underflow
            
            # Convert from log space
            probs[i] = np.exp(probs[i])
        
        probs /= np.max(probs)
        
        return probs
    
    def predict(self, X):
        probs = self.predict_prob(X)
        return self.classes[np.argmax(probs, axis=1)]
    
    def get_feature_distributions(self, X):
        """
        Get probability densities for both features across all classes
        
        Parameters:
        - feature_points: dict with 'similarity' and 'distance' arrays to evaluate densities
        
        Returns:
        - Dictionary with class-specific density values for both features
        """

        sim_X = np.linspace(X[:,0].min(), X[:,0].max(), 1000).reshape(-1, 1)
        dist_X = np.linspace(X[:,1].min(), X[:,1].max(), 1000).reshape(-1, 1)
        
        distributions = {
            'similarity': np.exp(self.similarity_kde.score_samples(sim_X)),
            'distance': np.exp(self.distance_kde.score_samples(dist_X))
        }
        
        return distributions


def nb_matrices(sim_matrix:np.ndarray, distance_matrix:np.ndarray, session_id):
    """
    Trains and fits a NB for a single pair of sessions.
    Parameters
    ----------
    sim_matrix : np.ndarray
        NxN matrix of similarity scores between units across two sessions.
    distance_matrix : np.ndarray
        NxN matrix of distance scores between units across two sessions.
    session_id : np.ndarray
        Vector of length N indicating the session ID for each unit (e.g. 0 for session 1, 1 for session 2).
    
    Returns
    -------
    np.ndarray
        NxN matrix of predicted probabilities for each unit pair being the same unit.
    """

    if sim_matrix.shape[0] != len(session_id) or distance_matrix.shape[0] != len(session_id):
        raise ValueError("Should have NxN matrices for similarity and distance, and a vector of length N for session IDs. N = number of units across both sessions.")
    
    N = len(session_id)
    labels = np.eye(N)
    within_session = (session_id[:, None] == session_id).astype(int)
    
    dist_within = distance_matrix[within_session==True].reshape(-1, 1)
    sim_within = sim_matrix[within_session==True].reshape(-1, 1)
    labels = labels[within_session==True].flatten()

    X_train = np.hstack((sim_within, dist_within))
    y_train = labels

    model = NaiveBayes(distance_bandwidth=1.0)
    model.fit(X_train, y_train)

    X_test = np.hstack((sim_matrix.reshape(-1, 1), distance_matrix.reshape(-1, 1)))
    predicted_probs = model.predict_prob(X_test)

    return predicted_probs.reshape(N, N)
