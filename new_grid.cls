netcdf grid_regular_SS {

    dimensions:

            grid_size = 12120 ;

            grid_corners = 4 ;

            grid_rank = 2 ;

    variables:

            int grid_dims(grid_rank) ;

            double grid_center_lat(grid_size) ;

                    grid_center_lat:units = "degrees" ;

                    grid_center_lat:bounds = "grid_corner_lat" ;

            double grid_center_lon(grid_size) ;

                    grid_center_lon:units = "degrees" ;

                    grid_center_lon:bounds = "grid_corner_lon" ;

            int grid_imask(grid_size) ;

                    grid_imask:units = "unitless" ;

                    grid_imask:coordinates = "grid_center_lon grid_center_lat" ;

            double grid_corner_lat(grid_size, grid_corners) ;

                    grid_corner_lat:units = "degrees" ;

            double grid_corner_lon(grid_size, grid_corners) ;

                    grid_corner_lon:units = "degrees" ;



    // global attributes:

                    :title = "Regular grid" ;

    }