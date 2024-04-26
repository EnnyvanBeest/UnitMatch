# import from plugins/custom_unitmatch_probability.py
"""Show how to add a custom unitmatch probability measure."""

from operator import itemgetter
import numpy as np

from phy import IPlugin
from phy.apps.template import from_sparse

def load_probability_templates(controller):
    try:
        out = controller.model._read_array(controller.model._find_path('probability_templates.npy'))
        out = np.atleast_2d(out)
        assert out.ndim == 2
        return out
    except IOError:
        return np.zeros((controller.n_templates, controller.n_templates))

class UnitMatchProbabilityPlugin(IPlugin):  
    def attach_to_controller(self, controller):
        print("Will set the similarity between clusters as UnitMatch's match probability!")

        # We cache this function in memory and on disk.
        # @controller.context.memcache ### not sure this is necessary, the fact that it's cached on disk makes it quite annoying if recomputing
        def um_probability(cluster_id):
            """Return the list of similar clusters to a given cluster."""
            # Templates of the cluster.
            temp_i = np.nonzero(controller.get_template_counts(cluster_id))[0]
            # The probability of the cluster to match with each template.
            probability_templates = load_probability_templates(controller)
            sims = np.max(probability_templates[temp_i, :], axis=0)

            def _sim_ij(cj):
                # Templates of the cluster.
                if cj < controller.model.n_templates:
                    return float(sims[cj])
                temp_j = np.nonzero(controller.get_template_counts(cj))[0]
                return float(np.max(sims[temp_j]))

            out = [(cj, _sim_ij(cj)) for cj in controller.supervisor.clustering.cluster_ids]
            return sorted(out, key=itemgetter(1), reverse=True)

        # We add the probability function.
        controller.similarity_functions['unitmatch_probability'] = um_probability

        # We set the probability function to the newly-defined one.
        controller.similarity = 'unitmatch_probability'