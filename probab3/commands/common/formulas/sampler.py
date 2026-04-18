import numpy as np
from probab3.commands.common.constants import *
from probab3.commands.common.general_code import LOG_FILE_NAME
from scipy import interpolate
from matplotlib import pyplot as plt
import logging

logging.basicConfig(filename=LOG_FILE_NAME, filemode='a', level=logging.INFO, format='[%(process)d]:[%(name)s]:%(levelname)s:%(message)s')
logger = logging.getLogger(__name__)

class Sampler():

    NUM_OF_POINTS_IN_DOMAIN = 20
    ZOOM_IN_THRESHOLD = 10

    def __init__(self, cdf, domains, *args):
        self.cdf = cdf
        self.domains = domains
        self.args = args

        self.domain_upper_boundary = domains[-1][-1]
        self.domain_lower_boundary = domains[0][0]

        self.total_cdf_volume = self.cdf(self.domain_upper_boundary, *args)
        self.cdf_points = []
        self.x_points = []

        if self.total_cdf_volume[0] == 0:
            raise CDFVolumeZero("The total cdf volume is zero. Cannot sample from this distribution.")
        
        last_cdf_value = 0
        num_of_empty_domains = 0 
        for domain in domains:
            threshold, num_of_points = self.get_threshold_num_of_points(len(domains) - num_of_empty_domains)
            cdf_points, x_points = self.get_invese_cdf_points_obey_threshold(domain, last_cdf_value, threshold, num_of_points)
            self.cdf_points += cdf_points
            self.x_points += list(x_points)
            if len(cdf_points) > 1:
                last_cdf_value = cdf_points[-1]
            if len(cdf_points) == 0:
                num_of_empty_domains += 1

            logger.debug(f"Got {len(cdf_points)} points for cdf construction in domain {domain}")

    def get_threshold_num_of_points(self, domains_num):
        threshold = np.ceil(self.ZOOM_IN_THRESHOLD / domains_num)
        threshold = threshold if threshold > 2 else self.ZOOM_IN_THRESHOLD
        num_of_points = np.ceil(self.NUM_OF_POINTS_IN_DOMAIN / domains_num)
        num_of_points = num_of_points if num_of_points > 4 else self.NUM_OF_POINTS_IN_DOMAIN
        return int(threshold), int(num_of_points)
    
    def get_invese_cdf_points_obey_threshold(self, domain, last_cdf_value=0, 
                                             threshold=ZOOM_IN_THRESHOLD, 
                                             num_of_points=NUM_OF_POINTS_IN_DOMAIN):
        
        cdf_points, xs = self.get_inverse_cdf_points(domain, last_cdf_value=last_cdf_value, num_of_points=num_of_points)
        if len(cdf_points) <= 1:
            # nothing of value in the domain
            logger.debug(f"Got no points from domain {domain} when constructing inverse cdf")
            return [], []
        
        while len(cdf_points) < threshold:
            logger.debug(f"Got {len(cdf_points)} cdf points < {threshold} zoom in threshold, in domain {domain}, zooming in")
            new_domain = (xs[0], xs[-1])
            cdf_points, xs = self.get_inverse_cdf_points(new_domain, last_cdf_value=last_cdf_value, num_of_points=num_of_points)
 
        return cdf_points, xs


    def get_inverse_cdf_points(self, domain, last_cdf_value=0, num_of_points=NUM_OF_POINTS_IN_DOMAIN):
        cdf_values = []

        x_values = np.linspace(domain[0], domain[-1], num_of_points)
        last_zero_index = 0
        first_one_index = len(x_values) - 1

        total_domain_cdf_value = self.cdf(x_values[-1], *self.args)[0]/self.total_cdf_volume[0]
        if total_domain_cdf_value == last_cdf_value:
            # nothing of value in the domain
            return [], []

        for index, x in enumerate(x_values):
            if index == 0:
               cdf_value = last_cdf_value 
            elif index == len(x_values) - 1:
                cdf_value = total_domain_cdf_value 
            else:
                try:
                    unnorm_cdf = self.cdf(x, *self.args)[0] 
                    cdf_value = unnorm_cdf/self.total_cdf_volume[0]
                except ZeroDivisionError as err:
                    logger.error(f"err={str(err)} at x={x}. Continuing")
                    cdf_values.append(cdf_values[-1] if cdf_values else last_cdf_value)
                    continue
                except Exception as err:
                    logger.error(f"err={str(err)} at x={x}")
                    raise err
                
            cdf_values.append(cdf_value)
            if cdf_value == 0:
                last_zero_index = index
            if cdf_value == 1:
                first_one_index = index
        
        return cdf_values[last_zero_index: first_one_index + 1], x_values[last_zero_index: first_one_index + 1]
    
    def evaluate_interp_inverse_cdf(self, points):
        logger.debug(f"cdf points {len(self.cdf_points)} x_points {len(self.x_points)}")
        interp_inv = interpolate.interp1d(self.cdf_points, self.x_points)
        return interp_inv(points)
    
    def plot_interp_inverse_cdf(self):
        p_range = np.linspace(0.00001, 1, 100)
        interp_inv = interpolate.interp1d(self.cdf_points, self.x_points, kind="slinear")
        interp_inv_cdf_values = interp_inv(p_range) 
        plt.plot(p_range, interp_inv_cdf_values)
        plt.scatter(self.cdf_points, self.x_points)
        plt.show()

    def plot_inverse_cdf(self):
        plt.scatter(self.cdf_points, self.x_points)
        plt.show()

    def sample(self, size=1):
        random_ps = np.random.uniform(size=size)
        return self.evaluate_interp_inverse_cdf(random_ps)

class MaxPDFValueNotFound(Exception):

    def __init__(self, message):
        super().__init__(message)

class CDFVolumeZero(Exception):

    def __init__(self, message):
        super().__init__(message)