import numpy as np
import pandas as pd
from scipy.optimize import minimize

class RadialAnalyzer:
    """
    Performs radial calculations on pre-processed spot and edge data.
    Processes the data on a per-file basis (grouped by File_ID), so that each
    imaging file (each circular colony) gets its own dynamic center.

    UPDATES:
      - Center calculation uses boundary-only fit (equal weight per angular sector)
      - Handles partial circles and asymmetric cell distributions
      - Diameter uses physical pattern sizes when known (offset_23/24/28/31)
      - Sign convention: NEGATIVE = toward center for radial_velocity, cos_theta
    """
    def __init__(self, spot_table: pd.DataFrame, edge_table: pd.DataFrame):
        """
        Parameters:
            spot_table (pd.DataFrame): Pre-processed DataFrame of spot data.
            edge_table (pd.DataFrame): Pre-processed DataFrame of edge data.
        """
        self.spot_table = spot_table
        self.edge_table = edge_table.copy()

    @staticmethod
    def _diameter_firstframe_hull(spots, time_col="SPOT_TIME", frame_col="FRAME", outer_pct=90.0):
        # pick earliest frame/time
        if time_col in spots.columns:
            t0 = np.nanmin(spots[time_col].values)
            sub = spots[np.isclose(spots[time_col].values, t0, rtol=0, atol=1e-9)]
        elif frame_col in spots.columns:
            f0 = int(np.nanmin(spots[frame_col].values))
            sub = spots[spots[frame_col] == f0]
        else:
            return np.nan, "firstframe_hull(no_time)"
        if sub.empty:
            return np.nan, "firstframe_hull(empty)"

        # rough center at t0
        cx0 = float(np.nanmedian(sub["POSITION_X"].values))
        cy0 = float(np.nanmedian(sub["POSITION_Y"].values))

        # outer rim points
        r = np.hypot(sub["POSITION_X"].values - cx0, sub["POSITION_Y"].values - cy0)
        if r.size < 3:
            return np.nan, "firstframe_hull(few)"
        rim = sub.loc[r >= np.percentile(r, outer_pct), ["POSITION_X","POSITION_Y"]].to_numpy()
        if rim.shape[0] < 3:
            return np.nan, "firstframe_hull(few_rim)"

        # hull width (fallback if SciPy missing)
        try:
            from scipy.spatial import ConvexHull, distance_matrix
            H = rim[ConvexHull(rim).vertices]
            D = float(distance_matrix(H, H).max())
        except Exception:
            P = rim
            dif = P[:,None,:] - P[None,:,:]
            D = float(np.sqrt((dif**2).sum(-1)).max())

        return D, "firstframe_hull"

    @staticmethod
    def boundary_only_fit(spots, n_sectors=36):
        """
        Fit circle using boundary-only approach with equal weight per angular sector.
        Works even with asymmetric cell distributions by sampling 1 point per direction.
        """
        x = spots['POSITION_X'].values
        y = spots['POSITION_Y'].values

        # Initial rough center
        cx0 = np.median(x)
        cy0 = np.median(y)

        # Divide into angular sectors
        angles = np.arctan2(y - cy0, x - cx0)
        sector_bins = np.linspace(-np.pi, np.pi, n_sectors + 1)

        # For each sector, find the OUTERMOST cell
        boundary_x = []
        boundary_y = []

        for i in range(n_sectors):
            sector_mask = (angles >= sector_bins[i]) & (angles < sector_bins[i+1])
            if np.any(sector_mask):
                sector_x = x[sector_mask]
                sector_y = y[sector_mask]
                # Get furthest cell from initial center in this sector
                distances = np.hypot(sector_x - cx0, sector_y - cy0)
                furthest_idx = np.argmax(distances)
                boundary_x.append(sector_x[furthest_idx])
                boundary_y.append(sector_y[furthest_idx])

        boundary_x = np.array(boundary_x)
        boundary_y = np.array(boundary_y)

        # Fit circle to boundary points
        def objective(params):
            cx, cy, radius = params
            distances = np.hypot(boundary_x - cx, boundary_y - cy)
            return np.var(distances)

        r0 = np.median(np.hypot(boundary_x - cx0, boundary_y - cy0))
        result = minimize(objective, [cx0, cy0, r0], method='Nelder-Mead')

        cx_fit, cy_fit, r_fit = result.x

        # Calculate fit quality
        distances = np.hypot(boundary_x - cx_fit, boundary_y - cy_fit)
        std_dev = np.std(distances)

        return float(cx_fit), float(cy_fit), float(r_fit * 2), std_dev  # Return diameter

    def process_file(self, file_id: str) -> pd.DataFrame:
        # --- slice tables ----------------------------------------------------
        spots = self.spot_table[self.spot_table['File_ID'] == file_id]
        edges = self.edge_table[self.edge_table['File_ID'] == file_id].copy()
        if spots.empty or edges.empty:
            print(f"Warning: No data for File_ID {file_id}")
            return edges

        # drop TrackMate's zero-length "edge 0→0"
        edges = edges[~((edges['SPOT_SOURCE_ID'] == edges['SPOT_TARGET_ID']) &
                        (edges['SPOT_SOURCE_ID'] == 0))].copy()

        # ====================================================================
        # CENTER: Boundary-only fit (handles asymmetric distributions)
        # ====================================================================
        center_x, center_y, diameter_boundary, fit_std = self.boundary_only_fit(spots)
        print(f"File {file_id}: Center (boundary fit): ({center_x:.2f}, {center_y:.2f}), " +
              f"Diameter: {diameter_boundary:.1f}µm (fit std: {fit_std:.1f})")

        # ====================================================================
        # DIAMETER: Use physical pattern size when known
        # ====================================================================
        PATTERN_SIZES = {
            'offset_23': 575,
            'offset_24': 700,
            'offset_28': 700,
            'offset_31': 575,
            'offset_26': 700,
            'offset_29': 700,
        }

        # Calculate backup diameters
        min_x = spots['POSITION_X'].min()
        max_x = spots['POSITION_X'].max()
        min_y = spots['POSITION_Y'].min()
        max_y = spots['POSITION_Y'].max()
        diameter_bbox_um = 0.5 * ((max_x - min_x) + (max_y - min_y))

        r = np.hypot(spots['POSITION_X'] - center_x, spots['POSITION_Y'] - center_y)
        r_clip = r[r <= np.percentile(r, 99.5)]
        diameter_p99_um = 2.0 * float(np.percentile(r_clip, 99)) if r_clip.size else np.nan

        # First-frame hull
        D_ff_um, src_ff = self._diameter_firstframe_hull(spots)

        # Decide on diameter source
        if 'offset_27' in file_id:
            # offset_27 has mixed pattern sizes - use boundary fit
            diameter_um = diameter_boundary
            diameter_source = "boundary_fit_mixed_pattern"
        else:
            # Check if we know physical pattern size
            physical_diameter = None
            for key, size in PATTERN_SIZES.items():
                if key in file_id:
                    physical_diameter = size
                    break

            if physical_diameter:
                diameter_um = physical_diameter
                diameter_source = "physical_pattern"
            else:
                # Use boundary fit as primary, with fallbacks
                if np.isfinite(diameter_boundary) and fit_std < 50:  # Good fit
                    diameter_um = diameter_boundary
                    diameter_source = "boundary_fit"
                elif np.isfinite(diameter_p99_um):
                    diameter_um = diameter_p99_um
                    diameter_source = "calculated_p99"
                else:
                    diameter_um = diameter_bbox_um
                    diameter_source = "calculated_bbox"

        if not np.isfinite(diameter_um) or diameter_um <= 0:
            diameter_um = diameter_bbox_um
            diameter_source = "fallback_bbox"

        print(f"  Final diameter: {diameter_um:.1f}µm (source: {diameter_source})")

        # --- collectors ------------------------------------------------------
        processed_idx = []
        u_vals, v_vals = [], []
        sx_vals, sy_vals, tx_vals, ty_vals = [], [], [], []
        r_vals, rv_vals, theta_vals = [], [], []
        cx_vals, cy_vals = [], []
        cos_vals = []
        angle_vals = []
        frame_interval_s = 600

        for idx, row in edges.iterrows():
            src = spots[spots['ID'] == row['SPOT_SOURCE_ID']]
            tgt = spots[spots['ID'] == row['SPOT_TARGET_ID']]
            if src.empty or tgt.empty:
                continue

            sx, sy = src[['POSITION_X', 'POSITION_Y']].values[0]
            tx, ty = tgt[['POSITION_X', 'POSITION_Y']].values[0]
            u, v   = tx - sx, ty - sy

            # --- time step (minutes) ---
            if 'SPOT_TIME' in spots.columns:
                ts_min = float(src['SPOT_TIME'].iloc[0]) * 60.0
                tt_min = float(tgt['SPOT_TIME'].iloc[0]) * 60.0
                dt = tt_min - ts_min
            elif 'POSITION_T' in spots.columns:
                ts_min = float(src['POSITION_T'].iloc[0]) / 60.0
                tt_min = float(tgt['POSITION_T'].iloc[0]) / 60.0
                dt = tt_min - ts_min
            else:
                dt = ((int(tgt['FRAME'].iloc[0]) - int(src['FRAME'].iloc[0])) * frame_interval_s) / 60.0

            if not np.isfinite(dt) or dt <= 0:
                continue

            # --- radial unit vector (inward) ---
            dx, dy = center_x - sx, center_y - sy
            r_source = np.hypot(dx, dy)
            if not np.isfinite(r_source) or r_source == 0:
                continue
            ux, uy = dx / r_source, dy / r_source

            # --- signed radial velocity (µm/min; NEGATIVE = toward center) ---
            rv = -(u * ux + v * uy) / dt

            eps = 1e-12
            mag = np.sqrt(u*u + v*v)
            if not np.isfinite(mag) or mag <= eps:
                continue

            # speed (µm/min)
            speed = mag / dt

            # directional bias (POSITIVE = toward center, same sign to match rose plots)
            cos_theta = -rv / speed

            # angle for rose plot: 0 = toward center, ±π = away from center
            vx, vy = u / mag, v / mag
            sin_inward = vx * uy - vy * ux
            cos_inward = vx * ux + vy * uy
            angle_to_center = np.arctan2(sin_inward, cos_inward)

            # collect
            processed_idx.append(idx)
            u_vals.append(u);               v_vals.append(v)
            sx_vals.append(sx);             sy_vals.append(sy)
            tx_vals.append(tx);             ty_vals.append(ty)
            r_vals.append(r_source);        rv_vals.append(rv)
            theta_vals.append(np.arctan2(uy, ux))
            cx_vals.append(center_x);       cy_vals.append(center_y)
            cos_vals.append(cos_theta)
            angle_vals.append(angle_to_center)

        # subset to processed rows BEFORE assigning vector columns
        edges = edges.loc[processed_idx].copy()

        edges['u']  = u_vals;                       edges['v']  = v_vals
        edges['source_x'] = sx_vals;                edges['source_y'] = sy_vals
        edges['target_x'] = tx_vals;                edges['target_y'] = ty_vals
        edges['radial_distance'] = r_vals
        edges['radial_velocity'] = rv_vals
        edges['radial_velocity_theta'] = theta_vals
        edges['cos_theta'] = cos_vals
        edges["angle_to_center"] = angle_vals
        edges['center_x'] = cx_vals;                edges['center_y'] = cy_vals

        # attach diameter scalars
        edges['diameter_um']            = diameter_um
        edges['diameter_source']        = diameter_source
        edges['diameter_bbox_um']       = diameter_bbox_um
        edges['diameter_pctl_um']       = diameter_p99_um
        edges['diameter_firstframe_um'] = D_ff_um
        edges['diameter_boundary_um']   = diameter_boundary

        return edges

    def process_all_files(self) -> pd.DataFrame:
        """
        Processes all files by grouping the data by 'File_ID', running the radial analysis
        for each file, and concatenating the results.

        Returns:
            A combined edge DataFrame with radial analysis columns for all files.
        """
        file_ids = self.spot_table['File_ID'].unique()
        processed_edges_list = []
        for fid in file_ids:
            processed_edges = self.process_file(fid)
            if not processed_edges.empty:
                processed_edges_list.append(processed_edges)
        if processed_edges_list:
            return pd.concat(processed_edges_list, ignore_index=True)
        else:
            return pd.DataFrame()

    def calculate_radial_persistence(self, edge_data: pd.DataFrame) -> pd.Series:
        edge_data['radial_persistence'] = edge_data['cos_theta']
        return edge_data['radial_persistence']
