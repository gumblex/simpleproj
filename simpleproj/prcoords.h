#ifndef PRCOORDS_H
#define PRCOORDS_H

PJ_LONLAT prcoords_wgs_gcj(PJ_LONLAT wgs, int check_china);
PJ_LONLAT prcoords_gcj_wgs(PJ_LONLAT gcj, int check_china);
PJ_LONLAT prcoords_gcj_bd(PJ_LONLAT gcj, int check_china);
PJ_LONLAT prcoords_bd_gcj(PJ_LONLAT bd, int check_china);
PJ_LONLAT prcoords_bd_wgs(PJ_LONLAT bd, int check_china);
PJ_LONLAT prcoords_wgs_bd(PJ_LONLAT bd, int check_china);
PJ_LONLAT prcoords_gcj_wgs_exact(PJ_LONLAT gcj, int check_china);
PJ_LONLAT prcoords_bd_gcj_exact(PJ_LONLAT bd, int check_china);
PJ_LONLAT prcoords_bd_wgs_exact(PJ_LONLAT bd, int check_china);

#endif // header
