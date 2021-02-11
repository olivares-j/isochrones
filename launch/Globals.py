#----------- Split ------------------------------
# Number of processes/parts to partition the data
size = 25
#------------------------------------------------

#--------------- Observables ----------------------------------
identifier   = "source_id"

gaia_mod     = ["parallax","BP","G","RP"]
gaia_obs     = ["parallax","bp","g","rp"]
gaia_unc     = ["parallax_error","bp_error","g_error","rp_error"]

dosmass_mod  = ["J","H","K"]
dosmass_obs  = ["Jmag","Hmag","Kmag"]
dosmass_unc  = ["e_Jmag","e_Hmag","e_Kmag"]

panstar_mod  = ["PS_g", "PS_r", "PS_i", "PS_z","PS_y"]
panstar_obs  = ["gmag","rmag","imag","zmag","ymag"]
panstar_unc  = ["e_gmag","e_rmag","e_imag","e_zmag","e_ymag"]

parameters   = ["age","mass","distance","AV"]
#---------------------------------------------------------------

#-------------- Transform lists -------------------------------
list_mod  = sum([gaia_mod,dosmass_mod,panstar_mod],[])
list_obs  = sum([gaia_obs,dosmass_obs,panstar_obs],[])
list_unc  = sum([gaia_unc,dosmass_unc,panstar_unc],[])

#-- Remove parallax from photometry -------
bands    = list_mod[1:]
phot_obs = list_obs[1:]
phot_unc = list_unc[1:]
#----------------------------------------------

bands_mag = [ band+"_mag" for band in bands]

columns_data = sum([[identifier],list_obs,list_unc],[])

stats_names = sum([[identifier],["mean_"+p for p in parameters],
					["std_"+p for p in parameters]],[])
#--------------------------------------------------------------

#------------ Files -------------------------------
dir_base   = "/raid/jromero/OCs/Perseus/" #Don't forget / at the end
dir_data   = dir_base + "data/"
dir_chain  = dir_base + "chains/"
dir_plots  = dir_base + "plots/"
dir_outs   = dir_base + "outputs/"
file_data  = dir_data + "members_test_7+2MASS+PanSTARRS.csv"
#---------------------------------------------------

#------- Miscellanea -------
nan_values = 99.0

# These are the filtering values of BP
label_BP = "bp"
limit_BP = 15.0

n_obs_min = 3   # Minimum number of observed bands
add_unc = 0.05  # Add uncertainty
nan_unc = 1.0   # If uncertainty is missing (but not the value)
#-------------------------------------------------

#-------- Av Prior ------------------------------------------
# Assumed Gaussian in Av units
prior_Av = {"loc":1.5,"scale":1.0,"lower":0.0,"upper":6.0}
#------------------------------------------------------------

#-------- Distance prior --------------------------------------
# Assumed Gaussian in parsec units
prior_distance = {"loc":310.,"scale":100,"lower":0,"upper":1000}
#--------------------------------------------------------------



