library(here)
source(here("upstream_scripts/01.4_ipr2lin.R"))
source(here("upstream_scripts/05a_ipr2da.R"))

## WP_020839904_full/WP_020839904
# ipr2lin
# ipr2lin(ipr = '../full_analysis_20210108/WP_020839904_full/WP_020839904.iprscan.tsv',
#         acc2info = '../full_analysis_20210108/WP_020839904_full/WP_020839904.acc2info.tsv',
#         suffix = '../full_analysis_20210108/WP_020839904_full/WP_020839904')
# # ipr2da, save as da, call next
# da <- ipr2da(infile_ipr = '../full_analysis_20210108/WP_020839904_full/WP_020839904.iprscan_cln.tsv',
#              prefix = '../full_analysis_20210108/WP_020839904_full/WP_020839904')
# # add ipr to blast file
# append_ipr(ipr_da = da, blast = '../full_analysis_20210108/WP_020839904_full/WP_020839904.cln.clust.tsv',
#            prefix = '../full_analysis_20210108/WP_020839904_full/WP_020839904')

# wp0208 <- read_tsv('../full_analysis_20210108/WP_020839904_full/WP_020839904.full_analysis.tsv') %>%
#  mutate(SLength = SLength.x) %>% select(-SLength.x, -SLength.y)
#
# write_tsv(wp0208, '../full_analysis_20210108/WP_020839904_full/WP_020839904.full_analysis.tsv')


## WP_043825948_full/WP_043825948
# ipr2lin
# ipr2lin(ipr = '../full_analysis_20210108/WP_043825948_full/WP_043825948.iprscan.tsv',
#         acc2info = '../full_analysis_20210108/WP_043825948_full/WP_043825948.acc2info.tsv',
#         suffix = '../full_analysis_20210108/WP_043825948_full/WP_043825948')
# # ipr2da, save as da, call next
# da <- ipr2da(infile_ipr = '../full_analysis_20210108/WP_043825948_full/WP_043825948.iprscan_cln.tsv',
#              prefix = '../full_analysis_20210108/WP_043825948_full/WP_043825948')
# # add ipr to blast file
# append_ipr(ipr_da = da, blast = '../full_analysis_20210108/WP_043825948_full/WP_043825948.cln.clust.tsv',
#            prefix = '../full_analysis_20210108/WP_043825948_full/WP_043825948')


## WP_096882215_full/WP_096882215
# ipr2lin
# ipr2lin(ipr = '../full_analysis_20210108/WP_096882215_full/WP_096882215.iprscan.tsv',
#         acc2info = '../full_analysis_20210108/WP_096882215_full/WP_096882215.acc2info.tsv',
#         suffix = '../full_analysis_20210108/WP_096882215_full/WP_096882215')
# # ipr2da, save as da, call next
# da <- ipr2da(infile_ipr = '../full_analysis_20210108/WP_096882215_full/WP_096882215.iprscan_cln.tsv',
#              prefix = '../full_analysis_20210108/WP_096882215_full/WP_096882215')
# add ipr to blast file
# append_ipr(ipr_da = da, blast = '../full_analysis_20210108/WP_096882215_full/WP_096882215.cln.clust.tsv',
#            prefix = '../full_analysis_20210108/WP_096882215_full/WP_096882215')

## WP_108717204_full/WP_108717204
# ipr2lin
# ipr2lin(ipr = '../full_analysis_20210108/WP_108717204_full/WP_108717204.iprscan.tsv',
#         acc2info = '../full_analysis_20210108/WP_108717204_full/WP_108717204.acc2info.tsv',
#         suffix = '../full_analysis_20210108/WP_108717204_full/WP_108717204')
# ipr2da, save as da, call next
# da1 <- ipr2da(infile_ipr = '../full_analysis_20210108/WP_108717204_full/WP_108717204.iprscan_cln.tsv',
#              prefix = '../full_analysis_20210108/WP_108717204_full/WP_108717204')
# # add ipr to blast file
# append_ipr(ipr_da = da1, blast = '../full_analysis_20210108/WP_108717204_full/WP_108717204.cln.clust.tsv',
#            prefix = '../full_analysis_20210108/WP_108717204_full/WP_108717204')

## WP_129996984_full/WP_129996984
# ipr2lin
# ipr2lin(ipr = '../full_analysis_20210108/WP_129996984_full/WP_129996984.iprscan.tsv',
#         acc2info = '../full_analysis_20210108/WP_129996984_full/WP_129996984.acc2info.tsv',
#         suffix = '../full_analysis_20210108/WP_129996984_full/WP_129996984')
# # ipr2da, save as da, call next
# da <- ipr2da(infile_ipr = '../full_analysis_20210108/WP_129996984_full/WP_129996984.iprscan_cln.tsv',
#              prefix = '../full_analysis_20210108/WP_129996984_full/WP_129996984')
# # add ipr to blast file
# append_ipr(ipr_da = da, blast = '../full_analysis_20210108/WP_129996984_full/WP_129996984.cln.clust.tsv',
#            prefix = '../full_analysis_20210108/WP_129996984_full/WP_129996984')

########## STAPH #############

## ABD20640_full/ABD20640
# # ipr2lin
# ipr2lin(ipr = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD20640_full/ABD20640.iprscan.tsv',
#         acc2info = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD20640_full/ABD20640.acc2info.tsv',
#         suffix = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD20640_full/ABD20640')
## ipr2da, save as da, call next
# da <- ipr2da(infile_ipr = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD20640_full/ABD20640.iprscan_cln.tsv',
#              prefix = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD20640_full/ABD20640')
# # add ipr to blast file
# append_ipr(ipr_da = da, blast = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD20640_full/ABD20640.cln.clust.tsv',
#            prefix = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD20640_full/ABD20640')

## ABD21022_full/ABD21022
# # ipr2lin
# lin <- ipr2lin(ipr = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD21022_full/ABD21022.iprscan.tsv',
#         acc2info = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD21022_full/ABD21022.acc2info.tsv',
#         suffix = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD21022_full/ABD21022')
# # ipr2da, save as da, call next
# da <- ipr2da(infile_ipr = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD21022_full/ABD21022.iprscan_cln.tsv',
#              prefix = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD21022_full/ABD21022')
# # add ipr to blast file
# fa <- append_ipr(ipr_da = da, blast = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD21022_full/ABD21022.cln.clust.tsv',
#            prefix = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD21022_full/ABD21022')

## ABD21741_full/ABD21741
# # ipr2lin
# lin <- ipr2lin(ipr = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD21741_full/ABD21741.iprscan.tsv',
#         acc2info = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD21741_full/ABD21741.acc2info.tsv',
#         suffix = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD21741_full/ABD21741')
# # ipr2da, save as da, call next
# da <- ipr2da(infile_ipr = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD21741_full/ABD21741.iprscan_cln.tsv',
#              prefix = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD21741_full/ABD21741')
# # add ipr to blast file
# fa <- append_ipr(ipr_da = da, blast = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD21741_full/ABD21741.cln.clust.tsv',
#            prefix = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD21741_full/ABD21741')

## ABD22038_full/ABD22038
# # ipr2lin
# lin <- ipr2lin(ipr = '../molevol_data/project_data/saureus/ABD22038_full/ABD22038_full/ABD22038.iprscan.tsv',
#         acc2info = '../molevol_data/project_data/saureus/ABD22038_full/ABD22038_full/ABD22038.acc2info.tsv',
#         suffix = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD22038_full/ABD22038')
# # ipr2da, save as da, call next
# da <- ipr2da(infile_ipr = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD22038_full/ABD22038.iprscan_cln.tsv',
#              prefix = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD22038_full/ABD22038')
# # add ipr to blast file
# fa <- append_ipr(ipr_da = da, blast = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD22038_full/ABD22038.cln.clust.tsv',
#                  prefix = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD22038_full/ABD22038')

## ABD22752_full/ABD22752
# # ipr2lin
# lin <- ipr2lin(ipr = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD22752_full/ABD22752.iprscan.tsv',
#         acc2info = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD22752_full/ABD22752.acc2info.tsv',
#         suffix = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD22752_full/ABD22752')
# # ipr2da, save as da, call next
# da <- ipr2da(infile_ipr = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD22752_full/ABD22752.iprscan_cln.tsv',
#              prefix = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD22752_full/ABD22752')
# # add ipr to blast file
# fa <- append_ipr(ipr_da = da, blast = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD22752_full/ABD22752.cln.clust.tsv',
#                  prefix = '../molevol_data/project_data/saureus/sausa300_0200-04/ABD22752_full/ABD22752')
#

lin <- ipr2lin(
  ipr = "../molevol_data/project_data/saureus/sausa300_0200-04/saureus_0204_quick_out/saureus_0204.iprscan.tsv",
  acc2info = "../molevol_data/project_data/saureus/sausa300_0200-04/saureus_0204_quick_out/saureus_0204.acc2info.tsv",
  suffix = "../molevol_data/project_data/saureus/sausa300_0200-04/saureus_0204_quick_out/saureus_0204"
)
