/*
 * Copyright (C) 2025 André Löfgren
 *
 * This file is part of Biceps.
 *
 * Biceps is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Biceps is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Biceps. If not, see <https://www.gnu.org/licenses/>.
 */
#pragma once

#define SEC_PER_YEAR 3.15576e7
#define PA_TO_MPA 1e-6

namespace MESH1D {
    enum DOMAIN_IDS {
        EMPTY_ID = 0,
        INTERIOR_ID = 1,
        WEST_ID = 2,
        EAST_ID = 4,
        BOUNDARY_ID = WEST_ID | EAST_ID,
        DOMAIN_ID = BOUNDARY_ID | INTERIOR_ID
    };
}

namespace MESH2D {
    enum CELL_TYPE {
        QUADRILATERAL,
        TRIANGLE_LEFT,
        TRIANGLE_RIGHT
    };

    enum DOMAIN_IDS {
        EMPTY_ID = 0,
        INTERIOR_ID = 1,
        NORTH_ID = 2,
        NORTH_WEST_ID = 4,
        WEST_ID = 8,
        SOUTH_WEST_ID = 16,
        SOUTH_ID = 32,
        SOUTH_EAST_ID = 64,
        EAST_ID = 128,
        NORTH_EAST_ID = 256,
        SURFACE_ID = NORTH_WEST_ID | NORTH_ID | NORTH_EAST_ID,
        BED_ID = SOUTH_WEST_ID | SOUTH_ID | SOUTH_EAST_ID,
        BOUNDARY_ID = NORTH_ID | NORTH_WEST_ID | WEST_ID | SOUTH_WEST_ID |
            SOUTH_ID | SOUTH_EAST_ID | EAST_ID | NORTH_EAST_ID,
        DOMAIN_ID = BOUNDARY_ID | INTERIOR_ID
    };
}

enum FSSA_VERSION {
    FSSA_NONE,
    FSSA_VERTICAL,
    FSSA_NORMAL 
};

enum TIMESTEP_SCHEME {
    EXPLICIT,
    SIMPLICIT,
    MIDPOINT,
    CRANK_NICOLSON
};

enum VELOCITY_COMPONENT {
    HORIZONTAL,
    VERTICAL
};

enum COORD_SYSTEM {
    CARTESIAN,
    NORMAL_TANGENTIAL
};
