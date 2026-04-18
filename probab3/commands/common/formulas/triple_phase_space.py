from probab3.commands.common.formulas.general import *
from probab3.commands.common.formulas.sampler import Sampler, MaxPDFValueNotFound, CDFVolumeZero
from probab3.commands.common.general_code import LOG_FILE_NAME
from scipy import integrate
import abc
import dataclasses
import scipy
import logging
from colorama import init as colorama_init
from colorama import Fore
from colorama import Style
from datetime import datetime
from probab3.commands.common.formulas.compiled_formulas import hyperbollic_dsigma_dEBdLBdCB, elliptic_dsigma_dEBdLBdCB


logging.basicConfig(filename=LOG_FILE_NAME, filemode='a', level=logging.INFO, format='[%(process)d]:[%(name)s]:%(levelname)s:%(message)s')
logger = logging.getLogger(__name__)

class ConeDist(abc.ABC):
    EB_SMALL_SHIFT = (10**(-7))
    CB_SMALL_SHIFT = (10**(-14))

    def __init__(self, ma, mb, ms, E0, L0, alpha):
        self.ma = ma
        self.mb = mb
        self.ms = ms
        self.E0 = E0
        self.L0 = L0
        self.alpha = alpha
        self.M = M(self.ma, self.mb, self.ms)
        self.mB = mB(self.ma, self.mb)
        self.m = m(ma, mb, ms)
        self.CB_sampler = None  
        self.over_all_int_value = -1
        self.max_pdf = None
        self.expected_EB = None
        self.expected_LB = None

    @abc.abstractmethod
    def dsigma_dEBdLBdCB(self, EB, LB, CB):
        pass

    @abc.abstractmethod
    def bounds_EB(self, CB):
        pass

    @abc.abstractmethod
    def general_bounds_EB(self):
        pass 

    def dsigma_dEBdLBdCB_normalized(self, EB, LB, CB):
        self.integrate_over_all() 
        return self.dsigma_dEBdLBdCB(EB, LB, CB)/self.over_all_int_value
    
    def EB_upper_bound(self):
        return (self.alpha*self.ma*self.mb*self.E0)/(self.ms*self.mB + self.alpha*self.ma*self.mb)
    
    def LB_upper_bound(self, EB):
        return (((G**2)*((self.ma*self.mb)**3))/(2*(-EB)*self.mB))**(1/2)

    def EB_verge_plus_inner(self, CB):
        return (G*(2*self.M*self.m + self.alpha*self.ma*self.mb) + 
                ((G*(2*self.M*self.m + self.alpha*self.ma*self.mb))**2 - 
                 (8*(self.L0**2)*((CB**2) - 1)*self.E0)/self.m)**(1/2)
                 )/((4*(self.L0**2)*((CB**2) - 1))/(self.m*self.alpha*self.ma*self.mb*G))

    def EB_verge_minus_inner(self, CB):
        return (G*(2*self.M*self.m + self.alpha*self.ma*self.mb) - 
                ((G*(2*self.M*self.m + self.alpha*self.ma*self.mb))**2 - 
                 (8*(self.L0**2)*(CB**2 - 1)*self.E0)/self.m)**(0.5)
                 )/((4*(self.L0**2)*(CB**2 - 1))/(self.m*self.alpha*self.ma*self.mb*G))
    
    def CB_verge_minus_inner_inner(self):
        return -np.sqrt(1 + (1/(8*(self.L0**2)*self.E0))*(G**2)*(
            (2*self.M*self.m + self.alpha*self.ma*self.mb)**2)*self.m, dtype=complex)

    def CB_verge_plus_inner_inner(self):
        return np.sqrt(1 + (1/(8*(self.L0**2)*self.E0))*(G**2)*(
            (2*self.M*self.m + self.alpha*self.ma*self.mb)**2)*self.m, dtype=complex)
    
    def LB_verge_minus(self, CB, EB):
        current_R = R(self.ma, self.mb, EB, self.alpha)
        return self.L0*CB - np.sqrt(((self.L0**2)*(CB**2 - 1) 
                              + 2*current_R*G*self.M*(self.m**2) 
                              + 2*(self.E0 - EB)*(current_R**2)*self.m), dtype=complex)

    def LB_verge_plus(self, CB, EB):
        current_R = R(self.ma, self.mb, EB, self.alpha)
        return self.L0*CB + np.sqrt(((self.L0**2)*(CB**2 - 1) 
                              + 2*current_R*G*self.M*(self.m**2) 
                              + 2*(self.E0 - EB)*(current_R**2)*self.m),
                              dtype=complex)

    def bounds_CB_negative(self):
        return [-1 + self.CB_SMALL_SHIFT, np.real(self.CB_verge_minus_inner_inner())]
        
    def bounds_CB_positive(self):
        return [np.real(self.CB_verge_plus_inner_inner()), 1]
        
    def bounds_LB(self, EB, CB):
        return [np.real(np.min([np.max([self.LB_verge_minus(CB, EB), 0]), self.LB_upper_bound(EB)])), 
                np.real(np.max([np.min([self.LB_upper_bound(EB), self.LB_verge_plus(CB, EB)]), 0]))]
    
    def integrate_over_all(self, override=False):
        if not override and self.over_all_int_value != -1:
            return self.over_all_int_value
        
        func = lambda x, y, z: self.dsigma_dEBdLBdCB(EB=y, LB=x, CB=z)
        
        negative_part = integrate.nquad(func, [self.bounds_LB, self.bounds_EB, self.bounds_CB_negative])
        positive_part = integrate.nquad(func, [self.bounds_LB, self.bounds_EB, self.bounds_CB_positive]) 
        self.over_all_int_value = negative_part[0] + positive_part[0]
        return self.over_all_int_value

    def get_EB_expected(self, override=False):
        if not override and self.expected_EB:
            return self.expected_EB
        
        over_all_int = self.integrate_over_all()
        
        func = lambda x, y, z: y*self.dsigma_dEBdLBdCB(EB=y, LB=x, CB=z)
        
        negative_part = integrate.nquad(func, [self.bounds_LB, self.bounds_EB, self.bounds_CB_negative])
        positive_part = integrate.nquad(func, [self.bounds_LB, self.bounds_EB, self.bounds_CB_positive]) 
        self.expected_EB = (negative_part[0] + positive_part[0])/over_all_int 
        return self.expected_EB

    def get_k(self, override=False):
        if not self.expected_EB or override:
            self.get_EB_expected(override=override)
        return self.expected_EB/self.E0

    def get_LB_expected(self, override=False):
        if not override and self.expected_LB:
            return self.expected_LB
        
        over_all_int = self.integrate_over_all()
        
        func = lambda x, y, z: x*self.dsigma_dEBdLBdCB(EB=y, LB=x, CB=z)
        
        negative_part = integrate.nquad(func, [self.bounds_LB, self.bounds_EB, self.bounds_CB_negative])
        positive_part = integrate.nquad(func, [self.bounds_LB, self.bounds_EB, self.bounds_CB_positive]) 
        self.expected_LB = (negative_part[0] + positive_part[0])/over_all_int 
        return self.expected_LB

    def get_l(self, override=False):
        if self.expected_LB is None or override:
            self.get_LB_expected(override=override)
        return self.expected_LB/self.L0
            
    
    def integrate_over_EBLB_maxCB(self, CB_max):
        func = lambda x, y, z: self.dsigma_dEBdLBdCB(EB=y, LB=x, CB=z)
        
        if CB_max < np.real(self.CB_verge_minus_inner_inner()):
            def bounds_CB_negative():
                return [-1, CB_max]
            
            return integrate.nquad(func, [self.bounds_LB, self.bounds_EB, bounds_CB_negative])
 
        def bounds_CB_negative():
                CB_bound = 0
                CB_minus_inner_inner = self.CB_verge_minus_inner_inner()
                if np.imag(CB_minus_inner_inner) == 0:
                    CB_bound = np.real(CB_minus_inner_inner)

                return [-1, CB_bound]        

        def bounds_CB_positive():
            CB_bound = 0
            CB_plus_inner_inner = self.CB_verge_plus_inner_inner()
            if np.imag(CB_plus_inner_inner) == 0:
                CB_bound = np.real(CB_plus_inner_inner)
            
            return [CB_bound, CB_max]

        if CB_max < np.real(self.CB_verge_plus_inner_inner()):
            return integrate.nquad(func, [self.bounds_LB, self.bounds_EB, bounds_CB_negative])

        return self.integrate_negative_positive_with_bounds(func, bounds_negative=[self.bounds_LB, self.bounds_EB, bounds_CB_negative], 
                                                            bounds_positive=[self.bounds_LB, self.bounds_EB, bounds_CB_positive]) 

    def integrate_negative_positive_with_bounds(self, func, bounds_negative, bounds_positive):
        negative_part = integrate.nquad(func, bounds_negative)
        positive_part = integrate.nquad(func, bounds_positive) 
        return (negative_part[0] + positive_part[0], np.mean([negative_part[1], positive_part[1]]))


    @abc.abstractmethod
    def bounds_EB_EBmax(self, EB_max, CB):
        pass

    def integrate_over_LB_maxEB(self, EB_max, CB):

        func = lambda x, y: self.dsigma_dEBdLBdCB(EB=y, LB=x, CB=CB)  
        
        def bounds_LB(EB):
            return [np.min([np.max([self.LB_verge_minus(CB, EB), 0]), self.LB_upper_bound(EB)]), 
                    np.max([np.min([self.LB_upper_bound(EB), self.LB_verge_plus(CB, EB)]), 0])]
        
        def bounds_EB():
            return self.bounds_EB_EBmax(EB_max, CB)
        
        return integrate.nquad(func, [bounds_LB, bounds_EB])   

    def integrate_over_maxLB(self, LB_max, EB, CB):
        func = lambda x: self.dsigma_dEBdLBdCB(EB=EB, LB=x, CB=CB)  

        def bounds_LB():
            return [np.min([np.max([self.LB_verge_minus(CB, EB), 0]), self.LB_upper_bound(EB)]), 
                    np.max([np.min([self.LB_upper_bound(EB), self.LB_verge_plus(CB, EB), LB_max]), 0])]
        
        return integrate.nquad(func, [bounds_LB])


    def initalize_CB_sampler(self):
        self.CB_sampler = Sampler(self.integrate_over_EBLB_maxCB, [[-1, np.real(self.CB_verge_minus_inner_inner())], [np.real(self.CB_verge_plus_inner_inner()), 1]])

    def initalize_EB_sampler(self, CB):
        self.EB_sampler = Sampler(self.integrate_over_LB_maxEB, 
                                  [self.bounds_EB(float(CB))], float(CB))

    def initalize_LB_sampler(self, EB, CB):
        self.LB_sampler = Sampler(self.integrate_over_maxLB, 
                                  [self.bounds_LB(EB, float(CB))], EB, float(CB)) 

    def max_dsigma_dEBdLBdCB_normalized(self, override=False):

        if not override and self.max_pdf:
            return self.max_pdf
        
        negative_dsigma = lambda x: -np.real(self.dsigma_dEBdLBdCB_normalized(EB=x[1], LB=x[2], CB=x[0]))

        CB_minus = np.real(self.CB_verge_minus_inner_inner())
        EB_bounds = self.bounds_EB(CB_minus - (self.CB_SMALL_SHIFT))
        LB_bounds = self.bounds_LB(EB_bounds[-1], CB_minus - (self.CB_SMALL_SHIFT))
        initial_guess = (1 - self.CB_SMALL_SHIFT, EB_bounds[1], LB_bounds[1])
        bounds_negative = [(-1 + self.CB_SMALL_SHIFT, CB_minus - self.CB_SMALL_SHIFT), (None, 0), (0, None)]
        bounds_positive = [(np.real(self.CB_verge_plus_inner_inner()) + self.CB_SMALL_SHIFT, 1 - self.CB_SMALL_SHIFT), (None, 0), (0, None)]
        
        def domain(x):
            EB_bounds = self.bounds_EB(CB=x[0])
            EB_constraint_upper = EB_bounds[1] - x[1]
            EB_constraint_lower = x[1] - EB_bounds[0]
            LB_bounds = self.bounds_LB(EB=x[1], CB=x[0])
            LB_constraint_upper = LB_bounds[1] - x[2]
            LB_constraint_lower = x[2] - LB_bounds[0]

            return np.array([EB_constraint_upper, EB_constraint_lower, LB_constraint_upper,LB_constraint_lower])
        
        constraints = [{'type': 'ineq', 'fun': lambda x: domain(x)}]

        try: 
            ans_negative = scipy.optimize.minimize(negative_dsigma, x0=initial_guess, method='SLSQP', bounds=bounds_negative, constraints=constraints)
            ans_positive = scipy.optimize.minimize(negative_dsigma, x0=initial_guess, method='SLSQP', bounds=bounds_positive, constraints=constraints)
        except TypeError as err:
            logger.warning(f"recieved {err}")
            logger.warning(f"E0={self.E0}, L0={self.L0}")
            raise err

        neg_ans = ans_negative.fun 
        if ans_negative.success == False: 
            logger.debug(f"Negative\n {ans_negative}")
            neg_ans = 0

        pos_ans = ans_positive.fun 
        if ans_positive.success == False: 
            pos_ans = 0 
            logger.debug(f"Positive \n {ans_positive}")

        self.max_pdf = np.max([-neg_ans, -pos_ans])
        if self.max_pdf == 0:
            raise MaxPDFValueNotFound(f"Could not find maximum value of pdf with m1={self.ma}, m2={self.mb}, m3={self.ms}, E0={self.E0}, L0={self.L0}")

        return self.max_pdf


    def max_dsigma_dEBdLB_normalized(self, CB_value):
        
        negative_dsigma = lambda x: -np.real(self.dsigma_dEBdLBdCB_normalized(EB=x[0], LB=x[1], CB=CB_value))

        EB_bounds = self.bounds_EB(CB_value)
        LB_bounds = self.bounds_LB(EB_bounds[-1], CB_value)
        initial_guess = (EB_bounds[1], LB_bounds[1])
        bounds_general = [(EB_bounds[0], EB_bounds[1]), (0, None)]
        
        def domain(x):
            LB_bounds = self.bounds_LB(EB=x[0], CB=CB_value)
            LB_constraint_upper = LB_bounds[1] - x[1]
            LB_constraint_lower = x[1] - LB_bounds[0]

            return np.array([LB_constraint_upper,LB_constraint_lower])
        
        constraints = [{'type': 'ineq', 'fun': lambda x: domain(x)}]

        try: 
            ans_general = scipy.optimize.minimize(negative_dsigma, x0=initial_guess, method='SLSQP', bounds=bounds_general, constraints=constraints)
        except TypeError as err:
            logger.warning(f"recieved {err}")
            logger.warning(f"E0={self.E0}, L0={self.L0}, CB={CB_value}")
            raise err

        ans = ans_general.fun 

        if ans_general.success == False:
            raise MaxPDFValueNotFound(f"Could not find maximum value of pdf with m1={self.ma}, m2={self.mb}, m3={self.ms}, E0={self.E0}, L0={self.L0}, CB={CB_value}")

        return -ans

    def max_dsigma_dLB_normalized(self, CB_value, EB_value):
        
        negative_dsigma = lambda x: -np.real(self.dsigma_dEBdLBdCB_normalized(EB=EB_value, LB=x[0], CB=CB_value))

        LB_bounds = self.bounds_LB(EB_value, CB_value)
        initial_guess = (LB_bounds[1])
        bounds_general = [(LB_bounds[0], LB_bounds[1])]

        try: 
            ans_general = scipy.optimize.minimize(negative_dsigma, x0=initial_guess, method='SLSQP', bounds=bounds_general)
        except TypeError as err:
            logger.warning(f"recieved {err}")
            logger.warning(f"E0={self.E0}, L0={self.L0}, CB={CB_value}, EB={EB_value}")
            raise err

        ans = ans_general.fun 

        if ans_general.success == False:
            raise MaxPDFValueNotFound(f"Could not find maximum value of pdf with m1={self.ma}, m2={self.mb}, m3={self.ms}, E0={self.E0}, L0={self.L0}, CB={CB_value}, EB={EB_value}")

        return -ans
    

    def sample(self, size=1, rejection=True):
        fusion = False
        if rejection:
            try:
                samples = self.sample_rejection(size=size)
            except MaxPDFValueNotFound as err:
                logger.warning(f"{Fore.RED}Got error: {err} \n Trying with inverse cdf sampling{Style.RESET_ALL}")
                fusion = True
        
        if fusion:
            samples = []
            resample_size = size
            while len(samples) < size:
                current_samples, resample_size = self.sample_fusion(size=resample_size)
                samples = samples + current_samples
    
        if not rejection and len(samples) < size:
            samples = []
            resample_size = size
            while len(samples) < size:
                current_samples, resample_size = self.sample_inverse_cdf(size=resample_size)
                samples = samples + current_samples

        return samples


    def sample_inverse_cdf(self, size=1):
        samples = []
        resample_size = 0

        logger.debug(f"Sampling CB")
        if self.CB_sampler == None:
            self.initalize_CB_sampler()

        CB_values = self.CB_sampler.sample(size=size)
        for CB_value in CB_values:
            EB_value = None
            LB_value = None
            try:
                logger.debug(f"Sampling EB")
                EB_sampler = Sampler(self.integrate_over_LB_maxEB, 
                                  [self.bounds_EB(float(CB_value))], float(CB_value))
                
                EB_value = EB_sampler.sample()[0]
                logger.debug(f"Sampling LB")
                LB_sampler = Sampler(self.integrate_over_maxLB, 
                                  [self.bounds_LB(EB_value, float(CB_value))], EB_value, float(CB_value))
                
                LB_value = LB_sampler.sample()[0]
                samples.append((CB_value, EB_value, LB_value))

            except CDFVolumeZero as err:
                logger.warning(f"{Fore.RED}Got CDFVolumeZero error: {err} when sampling CB_value={CB_value}, EB_value={EB_value}, LB_value={LB_value}{Style.RESET_ALL}") 
                logger.warning(f"{Fore.RED}Adding +1 to re-sample size {resample_size}.{Style.RESET_ALL}")
                resample_size += 1
                continue
            except Exception as err:
                logger.error(f"{Fore.RED}Got error: {err} when sampling CB_value={CB_value}, EB_value={EB_value}, LB_value={LB_value}{Style.RESET_ALL}") 
                raise err
        
        return samples, resample_size

    def sample_fusion(self, size=1):
        samples = []
        resample_size = 0

        logger.info(f"Sampling fusion")
        
        if self.CB_sampler == None:
            logger.debug(f"Initalizing CB sampler")
            self.initalize_CB_sampler()

        logger.debug(f"Inverse-cdf sampling CB")
        CB_values = self.CB_sampler.sample(size=size)
        for CB_value in CB_values:
            EB_value = None
            LB_value = None
            try:
                logger.debug(f"Rejection sampling EB, LB with CB={CB_value}")
                rejection_sample = self.sample_rejection_EB_LB(CB_value, size=1, padding_num=10)
                samples = samples + rejection_sample
                sampled = False if len(rejection_sample) < 1 else True 
            except Exception as err:
                logger.warning(f"{Fore.RED}Got error: {err} when rejection sampling EB, LB with CB_value={CB_value}, continuing with inverse cdf..{Style.RESET_ALL}")
                sampled = False

            try:
                if not sampled:
                    logger.debug(f"Inverse-cdf sampling EB")
                    EB_sampler = Sampler(self.integrate_over_LB_maxEB, 
                                      [self.bounds_EB(float(CB_value))], float(CB_value))
                
                    EB_value = EB_sampler.sample()[0]

                    try:
                        logger.debug(f"Rejection sampling LB with EB={EB_value}, CB={CB_value}")
                        rejection_sample = self.sample_rejection_LB(CB_value, EB_value, size=1, padding_num=10)
                        samples = samples + rejection_sample
                        LB_sampled = False if len(rejection_sample) < 1 else True 
                    except Exception as err:
                        logger.warning(f"{Fore.RED}Got error: {err} when rejection sampling LB with CB_value={CB_value}, EB={EB_value} continuing with inverse cdf..{Style.RESET_ALL}")
                        LB_sampled = False
                    
                    if not LB_sampled:
                        logger.debug(f"Inverse-cdf sampling LB")
                        LB_sampler = Sampler(self.integrate_over_maxLB, 
                                        [self.bounds_LB(EB_value, float(CB_value))], EB_value, float(CB_value))
                
                        LB_value = LB_sampler.sample()[0]
                        samples.append((CB_value, EB_value, LB_value))

            except CDFVolumeZero as err:
                logger.warning(f"{Fore.RED}Got CDFVolumeZero error: {err} when sampling CB_value={CB_value}, EB_value={EB_value}, LB_value={LB_value}{Style.RESET_ALL}") 
                logger.warning(f"{Fore.RED}Adding +1 to re-sample size {resample_size}.{Style.RESET_ALL}")
                resample_size += 1
                continue
            except Exception as err:
                logger.error(f"{Fore.RED}Got error: {err} when sampling CB_value={CB_value}, EB_value={EB_value}, LB_value={LB_value}{Style.RESET_ALL}") 
                raise err
        
        return samples, resample_size

    def sample_rejection(self, size=1, padding_num=10):
        samples = []
        sample_padding_length = size*2*padding_num

        if self.over_all_int_value == -1:
            self.integrate_over_all()

        max_pdf_value = self.max_dsigma_dEBdLBdCB_normalized()
        pdf_value_samples = list(np.random.uniform(low=0.0, high=max_pdf_value, size=sample_padding_length))
        
        CB_values_negative_domain = np.random.uniform(low=-1.0, high=np.real(self.CB_verge_minus_inner_inner()), size=size*padding_num)
        CB_values_negative_domain = np.random.uniform(low=np.real(self.CB_verge_plus_inner_inner()), high=1.0, size=size*padding_num)
        CB_values = list(CB_values_negative_domain) + list(CB_values_negative_domain)

        for sample_padding_num, CB_value in enumerate(CB_values):
            EB_value_bounds = self.bounds_EB(CB_value)
            if EB_value_bounds[0] == EB_value_bounds[1]:
                continue

            EB_value = np.random.uniform(low=EB_value_bounds[0], high=EB_value_bounds[1], size=1)[0]
            LB_value_bounds = self.bounds_LB(EB_value, CB_value)
            try:
                if LB_value_bounds[0] == LB_value_bounds[1]:
                    continue 
                LB_value = np.random.uniform(low=LB_value_bounds[0], high=LB_value_bounds[1], size=1)[0] 
            
            except OverflowError as err:
                logger.error(f"got error {err}")
                logger.error(f"LB_value_bounds = {LB_value_bounds}")
                logger.error(f"EB_value_bounds = {EB_value_bounds}, E0={self.E0}")
                logger.error(f"EB_value = {EB_value}, CB_value = {CB_value}, L0={self.L0}")
                raise err

            pdf_value = self.dsigma_dEBdLBdCB_normalized(EB_value, LB_value, CB_value)
            #click.echo(f"samples {pdf_value_samples[sample_padding_num]} and pdf_value{pdf_value}")
            if pdf_value_samples[sample_padding_num] <= pdf_value:
                samples.append((CB_value, EB_value, LB_value)) 

                if len(samples) == size:
                    return samples
        
        if len(samples) < size:
            logger.debug(f"padding_num = {padding_num} was not enough, generated {len(samples)}/{size} samples")
            logger.debug(f"zooming in")
            new_size = size - len(samples) 
            leftover_samples = self.sample_rejection(size=new_size, padding_num=2*padding_num)
            return samples + leftover_samples 
        
        return samples
        
    def sample_rejection_EB_LB(self, CB_value, size=1, padding_num=10):
        
        if padding_num > 400:
            logger.debug(f"recieved padding num {padding_num} > 400, returning in recursion.")
            return []

        samples = []
        sample_padding_length = size*2*padding_num

        max_pdf_value = self.max_dsigma_dEBdLB_normalized(CB_value)
        pdf_value_samples = list(np.random.uniform(low=0.0, high=max_pdf_value, size=sample_padding_length))

        EB_value_bounds = self.bounds_EB(CB_value)

        EB_values = np.random.uniform(low=EB_value_bounds[0], high=EB_value_bounds[1], size=size*padding_num)

        for sample_padding_num, EB_value in enumerate(EB_values):

            LB_value_bounds = self.bounds_LB(EB_value, CB_value)
            try:
                if LB_value_bounds[0] == LB_value_bounds[1]:
                    continue 
                LB_value = np.random.uniform(low=LB_value_bounds[0], high=LB_value_bounds[1], size=1)[0] 
            
            except OverflowError as err:
                logger.error(f"got error {err}")
                logger.error(f"LB_value_bounds = {LB_value_bounds}")
                logger.error(f"EB_value_bounds = {EB_value_bounds}, E0={self.E0}")
                logger.error(f"EB_value = {EB_value}, CB_value = {CB_value}, L0={self.L0}")
                raise err

            pdf_value = self.dsigma_dEBdLBdCB_normalized(EB_value, LB_value, CB_value)
            #click.echo(f"samples {pdf_value_samples[sample_padding_num]} and pdf_value{pdf_value}")
            if pdf_value_samples[sample_padding_num] <= pdf_value:
                samples.append((CB_value, EB_value, LB_value)) 

                if len(samples) == size:
                    return samples
        
        if len(samples) < size:
            logger.debug(f"padding_num = {padding_num} was not enough, generated {len(samples)}/{size} samples")
            logger.debug(f"zooming in")
            new_size = size - len(samples) 
            leftover_samples = self.sample_rejection_EB_LB(CB_value, size=new_size, padding_num=2*padding_num)
            return samples + leftover_samples 
        
        return samples


    def sample_rejection_LB(self, CB_value, EB_value, size=1, padding_num=10):
        if padding_num > 400:
            logger.warning(f"recieved padding num {padding_num} > 400, returning in recursion.")
            return []

        samples = []
        sample_padding_length = size*2*padding_num

        max_pdf_value = self.max_dsigma_dLB_normalized(CB_value, EB_value)
        pdf_value_samples = list(np.random.uniform(low=0.0, high=max_pdf_value, size=sample_padding_length))

        LB_value_bounds = self.bounds_LB(EB_value, CB_value)

        LB_values = np.random.uniform(low=LB_value_bounds[0], high=LB_value_bounds[1], size=size*padding_num)

        for sample_padding_num, LB_value in enumerate(LB_values):

            pdf_value = self.dsigma_dEBdLBdCB_normalized(EB_value, LB_value, CB_value)
            #click.echo(f"samples {pdf_value_samples[sample_padding_num]} and pdf_value{pdf_value}")
            if pdf_value_samples[sample_padding_num] <= pdf_value:
                samples.append((CB_value, EB_value, LB_value)) 

                if len(samples) == size:
                    return samples
        
        if len(samples) < size:
            logger.debug(f"padding_num = {padding_num} was not enough, generated {len(samples)}/{size} samples")
            logger.debug(f"zooming in")
            new_size = size - len(samples) 
            leftover_samples = self.sample_rejection_LB(size=new_size, padding_num=2*padding_num)
            return samples + leftover_samples 
        
        return samples


class EllipticDist(ConeDist):

    def dsigma_dEBdLBdCB_imported(self, EB, LB, CB):
        return elliptic_dsigma_dEBdLBdCB(EB/MSun, LB/MSun, CB, 
                                            self.ma/MSun, self.mb/MSun, self.ms/MSun, 
                                            self.alpha, self.E0/MSun, self.L0/MSun)  

    def dsigma_dEBdLBdCB(self, EB, LB, CB):
        #Ls = Ls_from_LBL0(LB, CB, self.L0)
        Ls = np.sqrt(LB**2 - 2*LB*self.L0*CB + self.L0**2)
        #current_R = R(self.ma, self.mb, EB, self.alpha) 
        current_R = (self.alpha*G*self.ma*self.mb)/(-2*EB)
        heaviside_term = -(Ls**2 -2*current_R*G*self.M*(self.m**2)
                          +2*self.m*(EB-self.E0)*(current_R**2))
        if heaviside_term <= 0:
            return 0

        const_prefactor = (2*(np.pi**4)*(G**2)*(self.M**(5/2))
                           *self.mB
                           )/((self.ma*self.mb*self.ms)**(3/2))
        term1_energies = (EB-self.E0)**(3/2) 
        term1 = LB/(Ls*(term1_energies)*((-EB)**(3/2)))

        inner_arccos = (1-((2*current_R*(EB - self.E0))/(G*self.mB*self.ms))
             )/(np.sqrt(1-((2*self.M*(EB-self.E0)*(Ls**2))/((G**2)*((self.mB*self.ms)**3))), dtype=complex))

        arccos_term = np.arccos(inner_arccos, dtype=complex)

        inner_radical_term = ((-2*self.M*(EB - self.E0))/((G**2)*(self.ms**3)*(self.mB**3))
                              *(Ls**2 -2*current_R*G*self.M*(self.m**2)
                              +2*self.m*(EB-self.E0)*(current_R**2)))
        radical_term = np.sqrt(inner_radical_term, dtype=complex)
        
        result = const_prefactor*term1*(arccos_term - radical_term)
        return result
    
    def bounds_EB(self, CB):
        try:
            EB_verge_plus = self.EB_verge_plus_inner(CB)
            EB_verge_minus = self.EB_verge_minus_inner(CB) 
            if np.imag(EB_verge_plus) and np.abs(np.imag(EB_verge_plus)) < np.abs(self.EB_SMALL_SHIFT*np.real(EB_verge_plus)):
                # neglect roundoff error imaginary parts 
                logger.warning(f"err in EB bounds: EB_verge_plus is complex {EB_verge_plus} but imaginary value is small, neglecting.")
                EB_verge_plus = np.real(EB_verge_plus) 
            if np.imag(EB_verge_minus) and np.abs(np.imag(EB_verge_minus)) < np.abs(self.EB_SMALL_SHIFT*np.real(EB_verge_minus)):
                # neglect roundoff error imaginary parts 
                logger.warning(f"err in EB bounds: EB_verge_minus is complex {EB_verge_minus} but imaginary value is small, neglecting.")
                EB_verge_minus = np.real(EB_verge_minus) 

            bounds = [np.max([self.E0, EB_verge_plus]), 
                np.min([self.EB_upper_bound(), EB_verge_minus])] 
        except Exception as err:
            logger.error(f"Got err in EB bounds: {err}")
            logger.error(f"CB = {CB}, E0 = {self.E0}, L0 = {self.L0}, EB_verge_plus {EB_verge_plus}, EB_verge_minus {EB_verge_minus}, EB_upper_bound {self.EB_upper_bound()}, CB_verge_minus_inner_inner {self.CB_verge_minus_inner_inner()}")
            raise err

        return bounds
    
    def bounds_EB_EBmax(self, EB_max, CB):
        return [np.max([self.E0, self.EB_verge_plus_inner(CB)]), 
                np.min([self.EB_upper_bound(), self.EB_verge_minus_inner(CB), EB_max])]

    def general_bounds_EB(self):
        return [self.E0, self.EB_upper_bound()]
        

    

class HyperbollicDist(ConeDist):

    def dsigma_dEBdLBdCB_imported(self, EB, LB, CB):
        return hyperbollic_dsigma_dEBdLBdCB(EB/MSun, LB/MSun, CB, 
                                            self.ma/MSun, self.mb/MSun, self.ms/MSun, 
                                            self.alpha, self.E0/MSun, self.L0/MSun) 
    
    def dsigma_dEBdLBdCB(self, EB, LB, CB):
        Ls = Ls_from_LBL0(LB, CB, self.L0)
        current_R = R(self.ma, self.mb, EB, self.alpha) 

        heaviside_term =(2*current_R*G*self.M*(self.m**2)
                          +2*self.m*(self.E0-EB)*(current_R**2)-Ls**2)
        
        if heaviside_term <= 0:
            return 0

        const_prefactor = (2*(np.pi**4)*(G**2)*(self.M**(5/2))
                           *self.mB
                           )/((self.ma*self.mb*self.ms)**(3/2))
        term1 = LB/(Ls*((self.E0-EB)**(3/2))*((-EB)**(3/2)))

        arccosh_inside = (1+((2*current_R*(self.E0-EB))/(G*self.mB*self.ms))
             )/((1+((2*self.M*(self.E0-EB)*(Ls**2))/((G**2)*((self.mB*self.ms)**3))))**(1/2)) 
        #arccosh_term = np.log(arccosh_inside +(((arccosh_inside**2) - 1)**(1/2)))
        arccosh_term = np.arccosh(arccosh_inside, dtype=complex)

        radical_term_inside = (((2*self.M*(self.E0-EB))/((G**2)*(self.ms**3)*(self.mB**3)))
                        *(2*current_R*G*self.M*(self.m**2)
                          +2*self.m*(self.E0-EB)*(current_R**2)-Ls**2))
        
        radical_term = np.sqrt(radical_term_inside, dtype=complex)
        
        return const_prefactor*term1*(radical_term - arccosh_term)

    def bounds_EB(self, CB):
        EB_verge_plus = self.EB_verge_plus_inner(CB)
        EB_verge_minus = self.EB_verge_minus_inner(CB)
        if np.imag(EB_verge_plus) and np.abs(np.imag(EB_verge_plus)) < np.abs(self.EB_SMALL_SHIFT*np.real(EB_verge_plus)):
                # neglect 10^-6 imaginary parts 
                logger.warning(f"err in EB bounds: EB_verge_plus is complex {EB_verge_plus} but imaginary value is small, neglecting.")
                EB_verge_plus = np.real(EB_verge_plus) 
        if np.imag(EB_verge_minus) and np.abs(np.imag(EB_verge_minus)) < np.abs(self.EB_SMALL_SHIFT*np.real(EB_verge_minus)):
                # neglect 10^-6 imaginary parts 
                logger.warning(f"err in EB bounds: EB_verge_minus is complex {EB_verge_minus} but imaginary value is small, neglecting.")
                EB_verge_minus = np.real(EB_verge_minus) 
        elif np.imag(EB_verge_plus) > 0:
            logger.error(f"complex EB bounds encountered.. EB verge plus: {EB_verge_plus}, using maximal bounds")
            EB_verge_plus = 10*self.E0 
            EB_verge_minus = (1 + self.EB_SMALL_SHIFT)*self.E0
         
        return [np.min([np.max([10*self.E0, EB_verge_plus]), (1 + self.EB_SMALL_SHIFT)*self.E0]), 
                np.min([(1 + self.EB_SMALL_SHIFT)*self.E0, EB_verge_minus])]

    def bounds_EB_EBmax(self, EB_max, CB):
        # TODO: not sure about what to do when this is complex..

        EB_verge_plus = self.EB_verge_plus_inner(CB)
        EB_verge_minus = self.EB_verge_minus_inner(CB)
        if np.imag(EB_verge_plus) and np.abs(np.imag(EB_verge_plus)) < np.abs(self.EB_SMALL_SHIFT*np.real(EB_verge_plus)):
                # neglect 10^-6 imaginary parts 
                logger.warning(f"err in EB bounds: EB_verge_plus is complex {EB_verge_plus} but imaginary value is small, neglecting.")
                EB_verge_plus = np.real(EB_verge_plus) 
        if np.imag(EB_verge_minus) and np.abs(np.imag(EB_verge_minus)) < np.abs(self.EB_SMALL_SHIFT*np.real(EB_verge_minus)):
                # neglect 10^-6 imaginary parts 
                logger.warning(f"err in EB bounds: EB_verge_minus is complex {EB_verge_minus} but imaginary value is small, neglecting.")
                EB_verge_minus = np.real(EB_verge_minus) 
        elif np.imag(EB_verge_plus) > 0:
            logger.error(f"complex EB bounds encountered.. EB verge plus: {EB_verge_plus}, using maximal bounds")
            EB_verge_plus = 10*self.E0 
            EB_verge_minus = (1 + self.EB_SMALL_SHIFT)*self.E0

        return [np.max([10*self.E0, EB_verge_plus]), 
                np.min([(1 + self.EB_SMALL_SHIFT)*self.E0, EB_verge_minus, EB_max])]

    def general_bounds_EB(self):
        return [10*self.E0, (1 + self.EB_SMALL_SHIFT)*self.E0]


class TripleSystem():

    def __init__(self, m1, m2, m3, E0, L0, alpha, calc_phase_space=True):
        self.E0 = E0
        self.L0 = L0
        self.alpha = alpha
        if m1 == m2 and m2 == m3:
            self.init_dist_equal(m1, m2, m3, E0, L0, alpha)
        else:
            self.init_dist_unequal(m1, m2, m3, E0, L0, alpha)

        self.P_dis = None

        self.hyperbollic_phase_space = None
        self.elliptic_phase_space = None
        self.fs_probabilities = None
        self.ims_probabilities = None

        if calc_phase_space:
            self.calculate_phase_spaces()

    def init_dist_unequal(self, m1, m2, m3, E0, L0, alpha):
        self.hyper123 = HyperbollicDist(m1, m2, m3, E0, L0, alpha)
        self.hyper231 = HyperbollicDist(m2, m3, m1, E0, L0, alpha)
        self.hyper312 = HyperbollicDist(m3, m1, m2, E0, L0, alpha)
        self.ellip123 = EllipticDist(m1, m2, m3, E0, L0, alpha)
        self.ellip231 = EllipticDist(m2, m3, m1, E0, L0, alpha)
        self.ellip312 = EllipticDist(m3, m1, m2, E0, L0, alpha)
        self.dists = [self.hyper123, self.hyper231, self.hyper312,
                      self.ellip123, self.ellip231, self.ellip312]
        
    def init_dist_equal(self, m1, m2, m3, E0, L0, alpha):
        self.hyper123 = HyperbollicDist(m1, m2, m3, E0, L0, alpha)
        self.hyper231 = self.hyper123
        self.hyper312 = self.hyper123
        self.ellip123 = EllipticDist(m1, m2, m3, E0, L0, alpha)
        self.ellip231 = self.ellip123
        self.ellip312 = self.ellip123
        self.dists = [self.hyper123, self.hyper231, self.hyper312,
                      self.ellip123, self.ellip231, self.ellip312]



    def calculate_phase_spaces(self):
        for dist in self.dists:
            logger.debug(f"{str(datetime.now())}::integrating over {dist}")
            dist.integrate_over_all() 
        
        self.hyperbollic_phase_space = self.hyper123.over_all_int_value + self.hyper231.over_all_int_value + self.hyper312.over_all_int_value
        self.elliptic_phase_space = self.ellip123.over_all_int_value + self.ellip231.over_all_int_value +self.ellip312.over_all_int_value
        self.P_dis = self.hyperbollic_phase_space / (self.hyperbollic_phase_space + self.elliptic_phase_space)


    def disintegration_probability(self):
        if not self.P_dis:
            self.P_dis = self.hyperbollic_phase_space / (self.hyperbollic_phase_space + self.elliptic_phase_space)
        return self.P_dis
    
    def m3_ims_probability(self):
        return self.ellip123.over_all_int_value / self.elliptic_phase_space

    def m3_fs_probability(self):
        return self.hyper123.over_all_int_value / self.hyperbollic_phase_space

    def get_fs_probabilities(self):
        if self.fs_probabilities is not None:
            return self.fs_probabilities
        dists = [self.hyper123, self.hyper231, self.hyper312]
        fs_probabilities = []
        for dist in dists:
            fs_probabilities.append(dist.over_all_int_value / self.hyperbollic_phase_space)

        self.fs_probabilities = fs_probabilities
        return fs_probabilities 

    def get_ims_probabilities(self):
        if self.ims_probabilities is not None:
            return self.ims_probabilities
        dists = [self.ellip123, self.ellip231, self.ellip312]
        ims_probabilities = []
        for dist in dists:
            ims_probabilities.append(dist.over_all_int_value / self.elliptic_phase_space)

        self.ims_probabilities = ims_probabilities
        return ims_probabilities

     
@dataclasses.dataclass            
class StateSample():
    ma: int
    mb: int
    CB: float
    EB: int
    LB: int

    def eB(self):
        return eB_from_EBLB(self.ma, self.mb, self.EB, self.LB)

