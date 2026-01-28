import os
import requests
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import rioxarray
from datetime import datetime
import pytz # Per gestire fusi orari

# --- CONFIGURAZIONE ---
OUTDIR = "output"
API_BASE = "https://radar-api.protezionecivile.it"
TZ_ROME = pytz.timezone('Europe/Rome') # Fuso orario Italiano

# Crea cartella output se non esiste
os.makedirs(OUTDIR, exist_ok=True)

# --- STILI E COLORI (INVARIATI) ---
colors_p = ["#ffffff","#bfe7f9","#7ed1f3","#00a6e6","#003f7b","#f4f89f","#e6ed3b",
            "#ffd800","#ff9500","#ff2f00","#b40a00","#840000","#dd007f"]
boundaries_p = [0,0.1,0.5,1,3,5,7,10,15,20,30,40,50]
cmap_p = ListedColormap(colors_p)
cmap_p.set_under("none")
norm_p = BoundaryNorm(boundaries_p, cmap_p.N, clip=False)

colors_p_cum = [
    "#ffffff", "#8fd3ff", "#00a6ff", "#0055ff", "#0000aa", 
    "#32cd32", "#008000", "#ffff00", "#ffcc00", "#ff9900", 
    "#ff4500", "#ff0000", "#b30000", "#ff1493", "#ff00ff", 
    "#9400d3", "#4b0082", "#dadada", "#909090", "#505050", "#000000"
]
boundaries_p_cum = [0,1,5,10,15,20,30,40,50,75,100,125,150,200,250,300,400,500,650,800,1000]
cmap_p_cum = ListedColormap(colors_p_cum)
cmap_p_cum.set_under("none")
norm_p_cum = BoundaryNorm(boundaries_p_cum, cmap_p_cum.N, clip=False)

# --- FUNZIONI API ---

def get_last_product_time(product_type):
    url = f"{API_BASE}/findLastProductByType?type={product_type}"
    try:
        r = requests.get(url, timeout=15)
        r.raise_for_status()
        data = r.json()
        if data.get('total', 0) > 0:
            return data['lastProducts'][0]['time']
    except Exception as e:
        print(f"[API] Errore check {product_type}: {e}")
    return None

def download_product(product_type, epoch_ms, filename):
    url_api = f"{API_BASE}/downloadProduct"
    payload = {"productType": product_type, "productDate": epoch_ms}
    try:
        r = requests.post(url_api, json=payload, timeout=15)
        r.raise_for_status()
        dl_url = r.json().get('url')
        if not dl_url: return False

        print(f"[DL] Scaricando {product_type}...")
        with requests.get(dl_url, stream=True, timeout=60) as r:
            r.raise_for_status()
            with open(filename, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        return True
    except Exception as e:
        print(f"[DL] Errore download {product_type}: {e}")
        return False

# --- UTILITY DATA ---

def get_formatted_filename(prefix, epoch_ms):
    """Genera nome file con timestamp locale: RADAR_SRI_20260128_1345.webp"""
    # Converti epoch UTC in datetime
    dt_utc = datetime.fromtimestamp(epoch_ms / 1000.0, pytz.utc)
    # Converti in locale (Roma)
    dt_local = dt_utc.astimezone(TZ_ROME)
    # Formatta stringa
    time_str = dt_local.strftime("%Y%m%d_%H%M")
    return f"{prefix}_{time_str}.webp", dt_local

# --- PLOTTING ---

def setup_map():
    fig = plt.figure(figsize=(10, 11))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([6, 19, 35, 47.5], crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE.with_scale('10m'), linewidth=1, zorder=5)
    ax.add_feature(cfeature.BORDERS.with_scale('10m'), linestyle=':', linewidth=1, zorder=5)
    gl = ax.gridlines(draw_labels=True, linestyle='--', alpha=0.3, zorder=6)
    gl.top_labels = False
    gl.right_labels = False
    return fig, ax

def save_map_webp(filename):
    plt.savefig(filename, format='webp', dpi=100, bbox_inches='tight', pil_kwargs={'quality': 70})
    plt.close()
    print(f"[OUT] Mappa salvata: {filename}")

def plot_sri_ir():
    print("--- Elaborazione SRI + IR ---")
    sri_time = get_last_product_time("SRI")
    if not sri_time: return

    # Nome file dinamico basato sul tempo del dato
    fname, dt_local = get_formatted_filename(f"{OUTDIR}/RADAR_SRI", sri_time)
    
    # Se il file esiste già, saltiamo (evita ri-elaborazione inutile)
    if os.path.exists(fname):
        print(f"[SKIP] File {fname} già esistente.")
        return

    ir_time = get_last_product_time("IR_108")
    f_sri, f_ir = "temp_sri.tif", "temp_ir.tif"
    
    if not download_product("SRI", sri_time, f_sri): return
    has_ir = False
    if ir_time: has_ir = download_product("IR_108", ir_time, f_ir)

    try:
        fig, ax = setup_map()

        if has_ir:
            da_ir = rioxarray.open_rasterio(f_ir, masked=True).squeeze()
            crs_ir = ccrs.Projection(da_ir.rio.crs)
            ax.pcolormesh(da_ir.x, da_ir.y, da_ir, transform=crs_ir, cmap='gray_r', alpha=0.6, zorder=1)

        da_sri = rioxarray.open_rasterio(f_sri, masked=True).squeeze()
        crs_sri = ccrs.Projection(da_sri.rio.crs)
        da_sri = da_sri.where(da_sri >= 0.1)

        cf = ax.pcolormesh(da_sri.x, da_sri.y, da_sri, transform=crs_sri, cmap=cmap_p, norm=norm_p, zorder=10)

        cbar = plt.colorbar(cf, ax=ax, orientation='horizontal', pad=0.05, aspect=30, extend='max')
        cbar.set_label("Precipitazione Istantanea (mm/h)")
        cbar.set_ticks(boundaries_p)

        # Titolo con ora locale
        plt.title(f"Radar SRI + Sat IR", loc='left', fontsize=12, fontweight='bold')
        plt.title(f"Valid: {dt_local.strftime('%d/%m/%Y %H:%M %Z')}", loc='right', fontsize=10)

        save_map_webp(fname)

    except Exception as e:
        print(f"Errore SRI: {e}")
    finally:
        if os.path.exists(f_sri): os.remove(f_sri)
        if os.path.exists(f_ir): os.remove(f_ir)

def plot_cum24():
    print("--- Elaborazione CUM 24H ---")
    cum_time = get_last_product_time("CUM24")
    if not cum_time: return

    fname, dt_local = get_formatted_filename(f"{OUTDIR}/RADAR_CUM24", cum_time)
    
    if os.path.exists(fname):
        print(f"[SKIP] File {fname} già esistente.")
        return

    f_cum = "temp_cum24.tif"
    if not download_product("CUM24", cum_time, f_cum): return

    try:
        fig, ax = setup_map()
        da_cum = rioxarray.open_rasterio(f_cum, masked=True).squeeze()
        da_cum = da_cum.where(da_cum >= 1)

        cf = ax.pcolormesh(da_cum.x, da_cum.y, da_cum, transform=ccrs.PlateCarree(),
                           cmap=cmap_p_cum, norm=norm_p_cum, zorder=10)

        cbar = plt.colorbar(cf, ax=ax, orientation='horizontal', pad=0.05, aspect=30, extend='max')
        cbar.set_label("Precipitazione 24h (mm)")
        cbar.set_ticks(boundaries_p_cum)

        plt.title(f"Precipitazione Cumulata 24h", loc='left', fontsize=12, fontweight='bold')
        plt.title(f"Valid: {dt_local.strftime('%d/%m/%Y %H:%M %Z')}", loc='right', fontsize=10)

        save_map_webp(fname)

    except Exception as e:
        print(f"Errore CUM24: {e}")
    finally:
        if os.path.exists(f_cum): os.remove(f_cum)

if __name__ == "__main__":
    plot_sri_ir()
    plot_cum24()
