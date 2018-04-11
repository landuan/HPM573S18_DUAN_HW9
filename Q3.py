import scr.RandomVariantGenerators as rndClasses
import scr.StatisticalClasses as Stat


TRANS_MATRIX = [
    [0.75,  0.15,   0,    0.1],   # WELL
    [0,     0,     1,    0],   # STROKE
    [0,     0.25,     0.55,  0.2],   # POST_STROKE
    ]

class HealthStats:
    """ health states of patients with risk of stroke """
    WELL = 0
    STROKE = 1
    P_STROKE = 2
    DEATH = 3

class Patient:
    def __init__(self, id):
        """ initiates a patient
        :param id: ID of the patient
        :param parameters: parameter object
        """

        self._id = id
        # random number generator for this patient
        self._rng = None
        self.healthstat=0
        self.survival=0

    def simulate(self, sim_length):
        """ simulate the patient over the specified simulation length """

        # random number generator for this patient
        self._rng = rndClasses.RNG(self._id)

        k = 0  # current time step

        # while the patient is alive and simulation length is not yet reached
        while self.healthstat!=3 and k  < sim_length:
            # find the transition probabilities of the future states
            trans_probs = TRANS_MATRIX[self.healthstat]
            # create an empirical distribution
            empirical_dist = rndClasses.Empirical(trans_probs)
            # sample from the empirical distribution to get a new state
            # (returns an integer from {0, 1, 2, ...})
            new_state_index = empirical_dist.sample(self._rng)
            # update health state
            self.healthstat =new_state_index[0]
            # increment time step
            k += 1
        self.survival=k+1

    def get_survival_time(self):
        """ returns the patient's survival time"""
        return self.survival


class Cohort():
    def __init__(self,id):
        self._initial_pop_size=2000
        self.survivaltime=[]
        self.id=id

    def simulate(self):
        for i in range(self._initial_pop_size):
            patient=Patient(self.id*self._initial_pop_size+i)
            patient.simulate(1000)
            self.survivaltime.append(patient.get_survival_time())

    def get_survival_time(self):
        return self.survivaltime


cohort=Cohort(1)
cohort.simulate()

sum_stat = Stat.SummaryStat("dsa",cohort.get_survival_time())
CI_of_Expected=sum_stat.get_t_CI(0.05)
Mean_Survival=sum_stat.get_mean()

print(Mean_Survival,CI_of_Expected)
print(cohort.get_survival_time())