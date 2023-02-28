
from pathlib import Path

from openalea.mtg import traversal, mtg
from openalea.plantgl.all import Scene
import math
from hydroshoot import architecture, display, model


def build_mtg(path_file: Path, is_show_scene: bool = True) -> (mtg.MTG, Scene):
    grapevine_mtg = architecture.vine_mtg(file_path=path_file)

    for v in traversal.iter_mtg2(grapevine_mtg, grapevine_mtg.root):
        architecture.vine_phyto_modular(grapevine_mtg, v)
        architecture.vine_mtg_properties(grapevine_mtg, v)
        architecture.vine_mtg_geometry(grapevine_mtg, v)
        architecture.vine_transform(grapevine_mtg, v)

    # Display of the plant mock-up (result in 'fig_01_plant_mock_up.png')
    mtg_scene = display.visu(grapevine_mtg, def_elmnt_color_dict=True, scene=Scene(), view_result=is_show_scene)
    return grapevine_mtg, mtg_scene

# %gui qt5

g, scene = build_mtg(Path('C:\GitModeles//virtualberry\VirtualBerry\HydroCall\grapevine_pot.csv'))
summary_results = model.run(g=g, wd=Path('C:\GitModeles//virtualberry\VirtualBerry\HydroCall'), scene=scene, gdd_since_budbreak=1000.)

# ================================================================================================

def partition(n,MTG=g):
    if ("in" or "cx") in g.node(n).label:
        return True
    else:
        return False


g.insert_scale(3,partition,"mm")

number=1

for lb in g.properties()['label']:
    if "mm" in g.properties()['label'][lb]:
        g.properties()['label'][lb]=g.properties()['label'][lb]+str(number)
        number+=1

# g.node(1).XX = g.node(1).baseXYZ[0]
# g.node(1).YY = g.node(1).baseXYZ[1]
# g.node(1).ZZ = g.node(1).baseXYZ[2]

for scale in range(1,g.max_scale()):
    print(scale)
    print(g.vertices(scale)[0])
    g.node(g.vertices(scale)[0]).XX = g.node(1).baseXYZ[0]
    g.node(g.vertices(scale)[0]).YY = g.node(1).baseXYZ[1]
    g.node(g.vertices(scale)[0]).ZZ = g.node(1).baseXYZ[2]


for element in g.vertices(4):
    g.node(element).XX = g.properties()['TopPosition'][element][0]
    g.node(element).YY = g.properties()['TopPosition'][element][1]
    g.node(element).ZZ = g.properties()['TopPosition'][element][2]
    if not g.node(g.parent(element)):
        g.node(element).length = math.sqrt((g.node(element).XX - g.node(1).XX)**2 +
                                           (g.node(element).YY - g.node(1).YY)**2 +
                                           (g.node(element).ZZ - g.node(1).ZZ)**2)
    else:
        g.node(element).length = math.sqrt((g.node(element).XX - g.node(g.parent(element)).XX) ** 2 +
                                           (g.node(element).YY - g.node(g.parent(element)).YY) ** 2 +
                                           (g.node(element).ZZ - g.node(g.parent(element)).ZZ) ** 2)


for element in g.vertices(4):
    g.node(element).C_demand = 1
    if len(g.components(g.complex(element))) == 3 and element == g.components(g.complex(element))[2]:
        g.node(element).C_supply = g.properties()['An'][element] * 12.01 *  g.properties()['leaf_area'][element]

for metamer in g.vertices(3):
    g.node(metamer).XX = g.properties()['XX'][g.components(metamer)[0]]
    g.node(metamer).YY = g.properties()['YY'][g.components(metamer)[0]]
    g.node(metamer).ZZ = g.properties()['ZZ'][g.components(metamer)[0]]
    g.node(metamer).TopDia = g.properties()['TopDiameter'][g.components(metamer)[0]]
    g.node(metamer).C_demand = 1
    if metamer == g.components(g.complex(metamer))[0]:
        g.node(g.complex(metamer)).XX = g.node(metamer).XX
        g.node(g.complex(metamer)).YY = g.node(metamer).YY
        g.node(g.complex(metamer)).ZZ = g.node(metamer).ZZ
    if len(g.components(metamer)) == 3:
        g.node(metamer).C_supply = g.properties()['An'][g.components(metamer)[2]]*12.01
    else:
        g.node(metamer).C_supply = 0
    if not g.node(g.parent(metamer)):
        g.node(metamer).length = math.sqrt((g.node(metamer).XX - g.node(1).XX)**2 +
                                           (g.node(metamer).YY - g.node(1).YY)**2 +
                                           (g.node(metamer).ZZ - g.node(1).ZZ)**2)
    else:
        g.node(metamer).length = math.sqrt((g.node(metamer).XX - g.node(g.parent(metamer)).XX) ** 2 +
                                           (g.node(metamer).YY - g.node(g.parent(metamer)).YY) ** 2 +
                                           (g.node(metamer).ZZ - g.node(g.parent(metamer)).ZZ) ** 2)


import Musca2022


Musca2022.bary_base(g,3)
Musca2022.vertices_semilength(g,3)

Dists=Musca2022.distances_baricentre(g,3)

Musca2022.C_allocation_SIMWAL2ses(g,1,3,Dists)

Dists=Musca2022.distances(g,3)