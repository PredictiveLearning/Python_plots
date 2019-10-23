

N_par=15
N_chains=360
Burn_In=0

plot_other_chains=1

chains_file='senna_g0'
chains_dir='./chains/'




fix_max_like=0
max_log_like=300.

plots_in_log=0

#Par_Names=['weight','Log_like',r'$\alpha_{\mathrm{SF}}$',r'$\alpha_{\mathrm{SF \_ burst}}$',r'$\beta_{\mathrm{burst}}$',
#           r'$k_{\mathrm{AGN}}$',r'$f_{\mathrm{BH}}$',r'$V_{\mathrm{BH}}$',
#           r'$\epsilon_{\mathrm{reheat}}$',r'$V_{\mathrm{reheat}}$',r'$\beta_{\mathrm{reheat}}$',
#           r'$\eta_{\mathrm{eject}}$',r'$V_{\mathrm{eject}}$',r'$\beta_{\mathrm{eject}}$',r'$\gamma_{\mathrm{reinc}}$',
 #          r'$Z_{\mathrm{hot}}$',r'$R_{\mathrm{merger}}$',r'$\alpha_{\mathrm{friction}}$',r'$M_{\mathrm{r.p.}}$']

Par_Names=['weight','Log_like',r'$\alpha_{\mathrm{SF}}$',r'$\alpha_{\mathrm{SF \_ burst}}$',r'$\beta_{\mathrm{burst}}$',
           r'$k_{\mathrm{AGN}}$',r'$f_{\mathrm{BH}}$',r'$V_{\mathrm{BH}}$',
           r'$\epsilon_{\mathrm{reheat}}$',r'$V_{\mathrm{reheat}}$',r'$\beta_{\mathrm{reheat}}$',
           r'$\eta_{\mathrm{eject}}$',r'$V_{\mathrm{eject}}$',r'$\beta_{\mathrm{eject}}$',r'$\gamma_{\mathrm{reinc}}$',
           r'$\alpha_{\mathrm{friction}}$',r'$M_{\mathrm{r.p.}}$']


one_one_size_small=[5,4]
one_two_size_small=[7,4]

one_two_size_large=[7,3]
one_three_size_large=[11,4]
one_four_size_large=[12,4]
one_five_size_large=[14,4]
one_three_size_large=[10,4]

two_one_size_small=[5,6]
two_two_size_small=[7,6]
two_two_size_large=[7,6]
two_four_size_large=[11,6]

three_one_size_small=[5,8]
three_two_size_small=[7,8]

Hen19_best_fit = [0.060, 0.50, 0.38, 0.0025, 0.066, 700, 5.6, 1.1e+02, 2.9, 5.5, 2.2e+02, 2, 1.2e10, 1.8, 5.1e+04]
Hen15_best_fit = [0.025, 0.60, 1.90, 0.0053, 0.041, 750, 2.6, 480, 0.72, 0.62, 100, 0.80, 3.0e10, 2.5, 1.2e+04]
Guo11_best_fit = [0.011, 0.56, 0.70, 0.0000, 0.030, 280, 4.0, 80., 3.2, 0.18, 90, 3.2, 0.0, 2.0, 0.0]

#MCMC_best_fit = [0.034, 0.11, 0.38, 0.00092, 0.066, 7e+02, 5.6, 1.1e+02, 2.9, 5.5, 2.2e+02, 2, 1.8e10, 1.8, 5.1e+04]
MCMC_best_fit = [0.073, 0.39, 1, 0.0044, 0.048, 5e+02, 6.4, 1.3e+02, 1.8, 8.7, 2.2e+02, 1.7, 1.4e+10, 2, 4.8e+04]