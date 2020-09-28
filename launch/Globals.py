#----------- Split ------------------------------
# Number of processes/parts to partition the data
size = 25
#------------------------------------------------

#--------------- Observables ----------------------------------
identifier   = "GDR2_ID"

gaia_mod     = ["parallax","BP","G","RP"]
gaia_obs     = ["parallax","BP","G","RP"]
gaia_unc     = ["parallax_error","e_BP","e_G","e_RP"]

dosmass_mod  = ["J","H","K"]
dosmass_obs  = ["J","H","Ks"]
dosmass_unc  = ["e_J","e_H","e_Ks"]

panstar_mod  = ["PS_g", "PS_r", "PS_i", "PS_z","PS_y"]
panstar_obs  = ["g_sdss","r_sdss","i_sdss","z_sdss","Y"]
panstar_unc  = ["e_g_sdss","e_r_sdss","e_i_sdss","e_z_sdss","e_Y"]

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
dir_base   = "/raid/jromero/OCs/USco/Isochrones"
dir_data   = dir_base + "/data/"
dir_chain  = dir_base + "/chains/"
dir_plots  = dir_base + "/plots/"
file_data  = dir_data + "members_GDR2_p05.csv"
#---------------------------------------------------

#------- Miscelanea -------
nan_values = 99.0

# These are the filtering values of BP
label_BP = "BP"
limit_BP = 15.0

n_obs_min = 3 # Minimum number of observed bands

add_unc = 0.05


