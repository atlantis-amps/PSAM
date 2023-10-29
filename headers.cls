netcdf SSM_v43_2018_0001 {
dimensions:
        nele = 25019 ;
        node = 16012 ;
        siglay = 10 ;
        siglev = 11 ;
        three = 3 ;
        time = UNLIMITED ; // (24 currently)
        DateStrLen = 26 ;
        maxnode = 12 ;
        maxelem = 10 ;
        four = 4 ;
        nobc = 87 ;
variables:
        int nprocs ;
                nprocs:long_name = "number of processors" ;
        int partition(nele) ;
                partition:long_name = "partition" ;
        float x(node) ;
                x:long_name = "nodal x-coordinate" ;
                x:units = "meters" ;
        float y(node) ;
                y:long_name = "nodal y-coordinate" ;
                y:units = "meters" ;
        float lon(node) ;
                lon:long_name = "nodal longitude" ;
                lon:standard_name = "longitude" ;
                lon:units = "degrees_east" ;
        float lat(node) ;
                lat:long_name = "nodal latitude" ;
                lat:standard_name = "latitude" ;
                lat:units = "degrees_north" ;
        float xc(nele) ;
                xc:long_name = "zonal x-coordinate" ;
                xc:units = "meters" ;
        float yc(nele) ;
                yc:long_name = "zonal y-coordinate" ;
                yc:units = "meters" ;
        float lonc(nele) ;
                lonc:long_name = "zonal longitude" ;
                lonc:standard_name = "longitude" ;
                lonc:units = "degrees_east" ;
        float latc(nele) ;
                latc:long_name = "zonal latitude" ;
                latc:standard_name = "latitude" ;
                latc:units = "degrees_north" ;
        float siglay(siglay, node) ;
                siglay:long_name = "Sigma Layers" ;
                siglay:standard_name = "ocean_sigma/general_coordinate" ;
                siglay:positive = "up" ;
                siglay:valid_min = -1.f ;
                siglay:valid_max = 0.f ;
                siglay:formula_terms = "sigma: siglay eta: zeta depth: h" ;
        float siglev(siglev, node) ;
                siglev:long_name = "Sigma Levels" ;
                siglev:standard_name = "ocean_sigma/general_coordinate" ;
                siglev:positive = "up" ;
                siglev:valid_min = -1.f ;
                siglev:valid_max = 0.f ;
                siglev:formula_terms = "sigma:siglay eta: zeta depth: h" ;
        float siglay_center(siglay, nele) ;
                siglay_center:long_name = "Sigma Layers" ;
                siglay_center:standard_name = "ocean_sigma/general_coordinate" ;
                siglay_center:positive = "up" ;
                siglay_center:valid_min = -1.f ;
                siglay_center:valid_max = 0.f ;
                siglay_center:formula_terms = "sigma: siglay_center eta: zeta_center depth: h_center" ;
        float siglev_center(siglev, nele) ;
                siglev_center:long_name = "Sigma Levels" ;
                siglev_center:standard_name = "ocean_sigma/general_coordinate" ;
                siglev_center:positive = "up" ;
                siglev_center:valid_min = -1.f ;
                siglev_center:valid_max = 0.f ;
                siglev_center:formula_terms = "sigma:siglay_center eta: zeta_center depth: h_center" ;
        float h_center(nele) ;
                h_center:long_name = "Bathymetry" ;
                h_center:standard_name = "sea_floor_depth_below_geoid" ;
                h_center:units = "m" ;
                h_center:positive = "down" ;
                h_center:grid = "grid1 grid3" ;
                h_center:coordinates = "latc lonc" ;
                h_center:grid_location = "center" ;
        float h(node) ;
                h:long_name = "Bathymetry" ;
                h:standard_name = "sea_floor_depth_below_geoid" ;
                h:units = "m" ;
                h:positive = "down" ;
                h:grid = "Bathymetry_Mesh" ;
                h:coordinates = "x y" ;
                h:type = "data" ;
        int nv(three, nele) ;
                nv:long_name = "nodes surrounding element" ;
        int iint(time) ;
                iint:long_name = "internal mode iteration number" ;
        float time(time) ;
                time:long_name = "time" ;
                time:units = "days since 1858-11-17 00:00:00" ;
                time:format = "modified julian day (MJD)" ;
                time:time_zone = "UTC" ;
        int Itime(time) ;
                Itime:units = "days since 1858-11-17 00:00:00" ;
                Itime:format = "modified julian day (MJD)" ;
                Itime:time_zone = "UTC" ;
        int Itime2(time) ;
                Itime2:units = "msec since 00:00:00" ;
                Itime2:time_zone = "UTC" ;
        char Times(time, DateStrLen) ;
                Times:time_zone = "PST" ;
        float zeta(time, node) ;
                zeta:long_name = "Water Surface Elevation" ;
                zeta:units = "meters" ;
                zeta:positive = "up" ;
                zeta:standard_name = "sea_surface_height_above_geoid" ;
                zeta:grid = "Bathymetry_Mesh" ;
                zeta:coordinates = "time lat lon" ;
                zeta:type = "data" ;
                zeta:location = "node" ;
        float dtfa(time, node) ;
                dtfa:long_name = "FSH + Water Height" ;
                dtfa:units = "meters" ;
                dtfa:positive = "down" ;
                dtfa:grid = "Bathymetry_Mesh" ;
                dtfa:coordinates = "time lat lon" ;
                dtfa:type = "data" ;
                dtfa:location = "node" ;
        float wts(time, siglev, node) ;
                wts:long_name = "Vertical Sigma Coordinate Velocity" ;
                wts:units = "s -1" ;
                wts:type = "data" ;
                wts:grid = "fvcom_grid" ;
        int nbe(three, nele) ;
                nbe:long_name = "elements surrounding each element" ;
        int ntsn(node) ;
                ntsn:long_name = "#nodes surrounding each node" ;
        int nbsn(maxnode, node) ;
                nbsn:long_name = "nodes surrounding each node" ;
        int ntve(node) ;
                ntve:long_name = "#elems surrounding each node" ;
        int nbve(maxelem, node) ;
                nbve:long_name = "elems surrounding each node" ;
        float a1u(four, nele) ;
                a1u:long_name = "a1u" ;
        float a2u(four, nele) ;
                a2u:long_name = "a2u" ;
        float aw0(three, nele) ;
                aw0:long_name = "aw0" ;
        float awx(three, nele) ;
                awx:long_name = "awx" ;
        float awy(three, nele) ;
                awy:long_name = "awy" ;
        float art2(node) ;
                art2:long_name = "Area of elements around a node" ;
        float art1(node) ;
                art1:long_name = "Area of Node-Base Control volume" ;
        float u(time, siglay, nele) ;
                u:long_name = "Eastward Water Velocity" ;
                u:standard_name = "eastward_sea_water_velocity" ;
                u:units = "meters s-1" ;
                u:grid = "fvcom_grid" ;
                u:type = "data" ;
                u:coordinates = "time siglay latc lonc" ;
                u:mesh = "fvcom_mesh" ;
                u:location = "face" ;
        float v(time, siglay, nele) ;
                v:long_name = "Northward Water Velocity" ;
                v:standard_name = "Northward_sea_water_velocity" ;
                v:units = "meters s-1" ;
                v:grid = "fvcom_grid" ;
                v:type = "data" ;
                v:coordinates = "time siglay latc lonc" ;
                v:mesh = "fvcom_mesh" ;
                v:location = "face" ;
        float tauc(time, nele) ;
                tauc:long_name = "bed stress magnitude from currents" ;
                tauc:note1 = "this stress is bottom boundary condtion on velocity field" ;
                tauc:note2 = "dimensions are stress/rho" ;
                tauc:units = "m^2 s^-2" ;
                tauc:grid = "fvcom_grid" ;
                tauc:type = "data" ;
                tauc:coordinates = "time latc lonc" ;
                tauc:mesh = "fvcom_mesh" ;
                tauc:location = "face" ;
        float uard_obcn(time, nobc) ;
                uard_obcn:long_name = "UARD at OBC" ;
                uard_obcn:units = "meters s-1" ;
                uard_obcn:grid = "fvcom_grid" ;
                uard_obcn:type = "data" ;
                uard_obcn:coordinates = "time obc" ;
                uard_obcn:mesh = "fvcom_mesh" ;
                uard_obcn:location = "obc" ;
        float xflux_obc(time, siglay, nobc) ;
                xflux_obc:long_name = "Xflux at OBC" ;
                xflux_obc:units = "tracer unit * m+3 s-1" ;
                xflux_obc:grid = "fvcom_grid" ;
                xflux_obc:type = "data" ;
                xflux_obc:coordinates = "time siglay obc" ;
                xflux_obc:mesh = "fvcom_mesh" ;
                xflux_obc:location = "obc" ;
        float ua(time, nele) ;
                ua:long_name = "Vertically Averaged x-velocity" ;
                ua:units = "meters s-1" ;
                ua:grid = "fvcom_grid" ;
                ua:type = "data" ;
        float va(time, nele) ;
                va:long_name = "Vertically Averaged y-velocity" ;
                va:units = "meters s-1" ;
                va:grid = "fvcom_grid" ;
                va:type = "data" ;
        float temp(time, siglay, node) ;
                temp:long_name = "temperature" ;
                temp:standard_name = "sea_water_temperature" ;
                temp:units = "degrees_C" ;
                temp:grid = "fvcom_grid" ;
                temp:coordinates = "time siglay lat lon" ;
                temp:type = "data" ;
                temp:mesh = "fvcom_mesh" ;
                temp:location = "node" ;
        float salinity(time, siglay, node) ;
                salinity:long_name = "salinity" ;
                salinity:standard_name = "sea_water_salinity" ;
                salinity:units = "1e-3" ;
                salinity:grid = "fvcom_grid" ;
                salinity:coordinates = "time siglay lat lon" ;
                salinity:type = "data" ;
                salinity:mesh = "fvcom_mesh" ;
                salinity:location = "node" ;
        float viscofm(time, siglay, nele) ;
                viscofm:long_name = "Horizontal Turbulent Eddy Viscosity For Momentum" ;
                viscofm:units = "m 2 s-1" ;
                viscofm:grid = "fvcom_grid" ;
                viscofm:coordinates = "x y" ;
                viscofm:type = "data" ;
        float viscofh(time, siglay, node) ;
                viscofh:long_name = "Horizontal Turbulent Eddy Viscosity For Scalars" ;
                viscofh:units = "m 2 s-1" ;
                viscofh:grid = "fvcom_grid" ;
                viscofh:coordinates = "x y" ;
                viscofh:type = "data" ;
        float km(time, siglev, node) ;
                km:long_name = "Turbulent Eddy Viscosity For Momentum" ;
                km:units = "m 2 s-1" ;
                km:grid = "fvcom_grid" ;
                km:coordinates = "x y" ;
                km:type = "data" ;
        float kh(time, siglev, node) ;
                kh:long_name = "Turbulent Eddy Viscosity For Scalars" ;
                kh:units = "m 2 s-1" ;
                kh:grid = "fvcom_grid" ;
                kh:coordinates = "x y" ;
                kh:type = "data" ;
        float kq(time, siglev, node) ;
                kq:long_name = "Turbulent Eddy Viscosity For Q2/Q2L" ;
                kq:units = "m 2 s-1" ;
                kq:grid = "fvcom_grid" ;
                kq:coordinates = "x y" ;
                kq:type = "data" ;
        float q2(time, siglev, node) ;
                q2:long_name = "Turbulent Kinetic Energy" ;
                q2:units = "m2 s-2" ;
                q2:grid = "fvcom_grid" ;
                q2:coordinates = "x y" ;
                q2:type = "data" ;
        float q2l(time, siglev, node) ;
                q2l:long_name = "Turbulent Kinetic Energy X Turbulent Macroscale" ;
                q2l:units = "m3 s-2" ;
                q2l:grid = "fvcom_grid" ;
                q2l:coordinates = "x y" ;
                q2l:type = "data" ;
        float l(time, siglev, node) ;
                l:long_name = "Turbulent Macroscale" ;
                l:units = "m3 s-2" ;
                l:grid = "fvcom_grid" ;
                l:coordinates = "x y" ;
                l:type = "data" ;
        int wet_nodes(time, node) ;
                wet_nodes:long_name = "Wet_Nodes" ;
                wet_nodes:grid = "fvcom_grid" ;
                wet_nodes:type = "data" ;
                wet_nodes:coordinates = "time lat lon" ;
                wet_nodes:mesh = "fvcom_mesh" ;
                wet_nodes:location = "node" ;
        int wet_cells(time, nele) ;
                wet_cells:long_name = "Wet_Cells" ;
                wet_cells:grid = "fvcom_grid" ;
                wet_cells:type = "data" ;
                wet_cells:coordinates = "time latc lonc" ;
                wet_cells:mesh = "fvcom_mesh" ;
                wet_cells:location = "face" ;
        int wet_nodes_prev_int(time, node) ;
                wet_nodes_prev_int:long_name = "Wet_Nodes_At_Previous_Internal_Step" ;
                wet_nodes_prev_int:grid = "fvcom_grid" ;
                wet_nodes_prev_int:type = "data" ;
                wet_nodes_prev_int:coordinates = "time lat lon" ;
                wet_nodes_prev_int:mesh = "fvcom_mesh" ;
                wet_nodes_prev_int:location = "node" ;
        int wet_cells_prev_int(time, nele) ;
                wet_cells_prev_int:long_name = "Wet_Cells_At_Previous_Internal_Step" ;
                wet_cells_prev_int:grid = "fvcom_grid" ;
                wet_cells_prev_int:type = "data" ;
                wet_cells_prev_int:coordinates = "time latc lonc" ;
                wet_cells_prev_int:mesh = "fvcom_mesh" ;
                wet_cells_prev_int:location = "face" ;
        int wet_cells_prev_ext(time, nele) ;
                wet_cells_prev_ext:long_name = "Wet_Cells_At_Previous_External_Step" ;
                wet_cells_prev_ext:grid = "fvcom_grid" ;
                wet_cells_prev_ext:type = "data" ;

// global attributes:
                :title = "SSM_2018_286_Riv_HRRR_dist_HSW_0.80_HCM_z0_0.6" ;
                :institution = "School for Marine Science and Technology" ;
                :source = "FVCOM_4.3" ;
                :history = "model started at: 05/10/2022   20:02" ;
                :references = "http://fvcom.smast.umassd.edu, http://codfish.smast.umassd.edu" ;
                :Conventions = "CF-1.0" ;
                :CoordinateSystem = "Cartesian" ;
                :CoordinateProjection = "proj=utm +ellps=WGS84 +zone=10" ;
                :Tidal_Forcing = "Tidal Forcing Time Series Title: SSM tidal elevation (NAVD88)" ;
                :River_Forcing = "THERE ARE 286 RIVERS IN THIS MODEL.\n",
                        "RIVER INFLOW IS ON THE nodes WHERE TEMPERATURE AND SALINITY ARE calculated IN THE MODEL.\n",
                        "THE FOLLOWING RIVER NAMES ARE USED:\n",
                        "Agate_East\n",
                        "Agate_West\n",
                        "Anderson_east\n",
                        "Anderson_west\n",
                        "Artondale\n",
                        "Burley_Cr\n",
                        "Butler_Cr\n",
                        "Campbell_Cr\n",
                        "Chambers_Cr\n",
                        "Coulter_Cr\n",
                        "Cranberry_Cr\n",
                        "Dana_Passage_North\n",
                        "Dana_Passage_South\n",
                        "Deer_Cr\n",
                        "Capitol_Lake\n",
                        "Capitol_Lake1\n",
                        "Dutcher_Cove\n",
                        "Ellis_Mission_Cr\n",
                        "Filucy_Bay\n",
                        "Fox_Island\n",
                        "Frye_Cove\n",
                        "Gallagher_Cove\n",
                        "Glen_Cove\n",
                        "Goldsborough_Cr\n",
                        "McCormick_Cr\n",
                        "Grant_East\n",
                        "Grant_West\n",
                        "Green_Cove\n",
                        "Gull_Harbor\n",
                        "Hale_Passage\n",
                        "Henderson_Inlet\n",
                        "Herron\n",
                        "Hope_Island\n",
                        "Jarrel_Cove\n",
                        "Johns_Cr\n",
                        "Kennedy_Schneider\n",
                        "Kennedy_Schneider1\n",
                        "Ketron\n",
                        "Ketron_Island\n",
                        "Mable_Taylor_Cr\n",
                        "Mayo_Cove\n",
                        "McAllister_Cr\n",
                        "McAllister_Cr1\n",
                        "McLane_Cove\n",
                        "McLane_Cr\n",
                        "McNeil_Isl\n",
                        "Mill_Cr\n",
                        "Minter_Cr\n",
                        "Moxlie_Cr\n",
                        "Nisqually_R\n",
                        "Nisqually_R1\n",
                        "Peale_Passage\n",
                        "Perry_Cr\n",
                        "Purdy_Cr\n",
                        "Rocky_Cr\n",
                        "Rosedale\n",
                        "Schneider_Cr\n",
                        "Sequalitchew_Cr\n",
                        "Sherwood_Cr\n",
                        "Skookum_Cr\n",
                        "Snodgrass_Cr\n",
                        "Squaxin_Island_East\n",
                        "Squaxin_Island_West\n",
                        "Sun_Pt\n",
                        "Tolmie\n",
                        "Van_Gelden\n",
                        "Vaughn\n",
                        "Whitman_Cove\n",
                        "Wilson_Pt\n",
                        "Woodard_Cr\n",
                        "Woodland_Cr\n",
                        "Young_Cove\n",
                        "Hylebos_Cr\n",
                        "Puyallup_R\n",
                        "Puyallup_R1\n",
                        "Bainbridge_Island_West\n",
                        "Blackjack_Cr\n",
                        "Chico_Cr\n",
                        "Dyes_Inlet\n",
                        "Gorst_Cr\n",
                        "Liberty_Bay\n",
                        "Miller_Bay\n",
                        "Green_R\n",
                        "Green_R1\n",
                        "Bainbridge_Island_East\n",
                        "Blake_Island\n",
                        "Buenna\n",
                        "Green_Valley_Cr\n",
                        "Curley_Cr\n",
                        "Des_Moines_Cr\n",
                        "Ellisport\n",
                        "Federal_Way\n",
                        "Gig_Harbor\n",
                        "Judd_Cr\n",
                        "Kitsap_NE\n",
                        "Lake_Washington\n",
                        "Lake_Washington1\n",
                        "Magnolia_Bch\n",
                        "Maury_Island\n",
                        "Miller_Cr\n",
                        "Olalla_Cr\n",
                        "Saltwater_St_Pk\n",
                        "Shingle_Mill_Cr\n",
                        "South_Snohomish\n",
                        "Tahlequah\n",
                        "University_Place\n",
                        "Dabob_Bay\n",
                        "Dosewallips_R\n",
                        "Dosewallips_R1\n",
                        "Duckabush_R\n",
                        "Duckabush_R1\n",
                        "Hamma_Hamma_R\n",
                        "Hamma_Hamma_R1\n",
                        "Kitsap_Hood\n",
                        "Lynch_Cove\n",
                        "Lynch_Cove1\n",
                        "NW_Hood\n",
                        "Port_Gamble\n",
                        "Quilcene\n",
                        "Quilcene1\n",
                        "Skokomish_R\n",
                        "Skokomish_R1\n",
                        "Tahuya\n",
                        "Tahuya1\n",
                        "Skagit_R\n",
                        "Skagit_R1\n",
                        "Snohomish_R\n",
                        "Snohomish_R1\n",
                        "Stillaguamish_R\n",
                        "Stillaguamish_R1\n",
                        "Whidbey_east\n",
                        "Port_Townsend\n",
                        "Whidbey_west\n",
                        "Birch_Bay\n",
                        "Lopez_Island\n",
                        "Nooksack_R\n",
                        "Nooksack_R1\n",
                        "Orcas_Island\n",
                        "Samish_Bell_south\n",
                        "Samish_Bell_south1\n",
                        "San_Juan_Island\n",
                        "Whatcom_Bell_north\n",
                        "Clallam_Bay\n",
                        "Discovery_Bay\n",
                        "Dungeness_R\n",
                        "Dungeness_R1\n",
                        "Elwha_R\n",
                        "Elwha_R1\n",
                        "North_Olympic\n",
                        "Port_Angeles\n",
                        "Sequim_Bay\n",
                        "Fraser_R\n",
                        "Fraser_R1\n",
                        "Howe_Sound\n",
                        "Sunshine_Coast\n",
                        "Vancouver_Isl_C\n",
                        "Vancouver_Isl_N\n",
                        "Vancouver_Isl_S\n",
                        "Victoria_SJdF\n",
                        "Boston_Harbor\n",
                        "Carlyon\n",
                        "Chambers_Creek\n",
                        "Fort_Lewis\n",
                        "Hartstene\n",
                        "LOTT\n",
                        "McNeil_Is\n",
                        "Rustlewood\n",
                        "Seashore_Villa\n",
                        "Shelton\n",
                        "Tamoshan\n",
                        "Taylor_Bay\n",
                        "Puyallup\n",
                        "Tacoma_Central\n",
                        "Tacoma_North\n",
                        "US_Oil_&_Refining\n",
                        "West_Rock\n",
                        "Bainbridge_Kitsap_Co_7\n",
                        "Bremerton\n",
                        "Central_Kitsap\n",
                        "Port_Orchard\n",
                        "Suquamish\n",
                        "Alderwood\n",
                        "Bainbridge_Island_City\n",
                        "Brightwater\n",
                        "Edmonds\n",
                        "Gig_Harbor\n",
                        "Kitsap_Co_Kingston\n",
                        "Lakota\n",
                        "Lynnwood\n",
                        "Manchester\n",
                        "Messenger_House\n",
                        "Midway\n",
                        "Miller_Creek\n",
                        "Redondo\n",
                        "Salmon_Creek\n",
                        "South_King\n",
                        "Vashon\n",
                        "West_Point\n",
                        "Alderbrook\n",
                        "Port_Gamble\n",
                        "Coupeville\n",
                        "Everett_Snohomish\n",
                        "Kimberly_Clark\n",
                        "La_Conner\n",
                        "Lake_Stevens_001\n",
                        "Lake_Stevens_002\n",
                        "Langley\n",
                        "Marysville\n",
                        "Mt_Vernon\n",
                        "Mukilteo\n",
                        "Oak_Harbor_Lagoon\n",
                        "Oak_Harbor_RBC\n",
                        "OF100\n",
                        "Penn_Cove\n",
                        "Skagit_County_2_Big_Lake\n",
                        "Snohomish\n",
                        "Stanwood\n",
                        "Swinomish\n",
                        "Tulalip\n",
                        "Warm_Beach_Campground\n",
                        "Port_Ludlow\n",
                        "Port_Townsend\n",
                        "Port_Townsend_Paper\n",
                        "Anacortes\n",
                        "Bellingham\n",
                        "Birch_Bay\n",
                        "Blaine\n",
                        "BP_Cherry_Point\n",
                        "Conoco_Phillips\n",
                        "Eastsound_Orcas_Village\n",
                        "Eastsound_Water_District\n",
                        "Fisherman_Bay\n",
                        "Friday_Harbor\n",
                        "Intalco\n",
                        "Larrabee_State_Park\n",
                        "Lummi_Goose_Pt\n",
                        "Lummi_Sandy_Pt\n",
                        "Makah\n",
                        "Roche_Harbor\n",
                        "Rosario_Utilities\n",
                        "Shell_Oil\n",
                        "Tesoro\n",
                        "Whidbey_Naval_Station\n",
                        "Clallam_Bay_POTW\n",
                        "Clallam_DOC\n",
                        "Nippon_Paper\n",
                        "Port_Angeles\n",
                        "Sekiu\n",
                        "Sequim\n",
                        "Annacis\n",
                        "Iona\n",
                        "Lions_Gate\n",
                        "Lulu\n",
                        "NW_Langley\n",
                        "Clover_Point\n",
                        "Gulf_Islands\n",
                        "Macaulay\n",
                        "Saanich\n",
                        "Cushman_No_2\n",
                        "Chehalis_R\n",
                        "Chehalis_R1\n",
                        "Willapa_R\n",
                        "Willapa_R1\n",
                        "Columbia_R\n",
                        "Columbia_R1\n",
                        "Willamette_R\n",
                        "Willamette_R1\n",
                        "Holberg\n",
                        "North_East_Vancouver_Island\n",
                        "Brooks_Peninsula\n",
                        "Tahsis\n",
                        "Nimpkish_River\n",
                        "Tsitika_River\n",
                        "Gold_River\n",
                        "Salmon_River\n",
                        "Clayoquot\n",
                        "Campbell_River\n",
                        "Quadra_Islands\n",
                        "Alberni_Inlet\n",
                        "Owikeno_Lake\n",
                        "Homathco_River\n",
                        "Klinaklini_River\n",
                        "Seymour_Inlet\n",
                        "Knight_Inlet\n",
                        "Toba_Inlet\n",
                        "Neil_Creek" ;
                :GroundWater_Forcing = "GROUND WATER FORCING IS OFF!" ;
                :Surface_Heat_Forcing = "FVCOM variable surface heat forcing file:\n",
                        "FILE NAME:ssm_HSHORTDist_dist_hfx_2018.nc\n",
                        "SOURCE:FVCOM grid (unstructured) surface forcing\n",
                        "Unknown start date meta data format" ;
                :Surface_Wind_Forcing = "FVCOM variable surface Wind forcing:\n",
                        "FILE NAME:ssm_HSHORTDist_dist_hfx_2018.nc\n",
                        "SOURCE:FVCOM grid (unstructured) surface forcing\n",
                        "Unknown start date meta data format" ;
                :Surface_PrecipEvap_Forcing = "SURFACE PRECIPITATION FORCING IS OFF" ;
}
