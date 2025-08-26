#!/usr/bin/env python3
"""
plip_xml_to_json.py

Usage:
    python plip_xml_to_json.py report.xml report.json
"""
import sys
import json
import xml.etree.ElementTree as ET
import re
from collections import defaultdict

def tag_localname(tag):
    """Return local tag name (strip namespace) lowercase."""
    return tag.split('}')[-1].lower()

def find_residue_info(elem):
    """
    Try to extract residue info dict {'resname', 'resnr', 'chain'} from an element
    Candidate sources:
      - element attributes like resname,resnr,chain,name,number
      - child elements <resname>, <resnr>, <residueName>, <residueNumber>, <chain>
      - text inside a child <residue> with attributes
    Returns dict or None if nothing found.
    """
    attrib = elem.attrib
    # try many attribute names
    name = (attrib.get('resname') or attrib.get('residue_name') or attrib.get('residue') or
            attrib.get('name') or attrib.get('resName'))
    num = (attrib.get('resnr') or attrib.get('residue_number') or attrib.get('residue_id') or
           attrib.get('number') or attrib.get('seq') or attrib.get('resSeq'))
    chain = (attrib.get('chain') or attrib.get('chain_id') or attrib.get('chainId'))

    # try child nodes if attributes not present
    if name is None or num is None:
        for child in elem:
            t = tag_localname(child.tag)
            text = child.text.strip() if child.text else None
            if t in ('resname', 'residue_name', 'residue', 'name', 'residueName') and text:
                name = name or text
            if t in ('resnr', 'residue_number', 'residue_id', 'number', 'residueNumber') and text:
                num = num or text
            if t in ('chain', 'chain_id', 'chainid') and text:
                chain = chain or text

    # if child <residue> contains attributes or children, recursively check those
    if (name is None or num is None):
        for child in elem:
            if 'resid' in tag_localname(child.tag) or 'resid' in (child.attrib.get('type','').lower() if 'type' in child.attrib else ''):
                # try attributes of child
                c_at = child.attrib
                name = name or c_at.get('resname') or c_at.get('name')
                num = num or c_at.get('resnr') or c_at.get('number')
                chain = chain or c_at.get('chain')

    if not name and not num:
        return None

    # normalize numeric residue number to int if possible
    try:
        rn = int(re.search(r"-?\d+", str(num)).group())
    except Exception:
        rn = str(num)

    out = {'resname': str(name) if name is not None else None,
           'resnr': rn,
           'chain': str(chain) if chain is not None else ""}
    return out

def walk_and_collect(root):
    """
    Walk the XML and return a dict with keys:
      'hydrogen_bonds', 'hydrophobic_contacts', 'ionic_interactions', 'water_bridges'
    Each mapped to a list of residue dicts as returned by find_residue_info.
    """
    mapping = {
        'hydrogen': 'hydrogen_bonds',
        'hbond': 'hydrogen_bonds',
        'hydrophobic': 'hydrophobic_contacts',
        'hydro': 'hydrophobic_contacts',
        'ionic': 'ionic_interactions',
        'salt': 'ionic_interactions',
        'water': 'water_bridges',
        'bridge': 'water_bridges'
    }

    out = defaultdict(list)

    for elem in root.iter():
        t = tag_localname(elem.tag)
        # check tag or type attribute for keywords
        candidate = None
        # first check tag
        for kw, outkey in mapping.items():
            if kw in t:
                candidate = outkey
                break
        # if not found, check 'type' or 'interactionType' attributes
        if candidate is None:
            typ = (elem.attrib.get('type') or elem.attrib.get('interaction') or elem.attrib.get('interactionType') or "")
            typ = typ.lower()
            for kw, outkey in mapping.items():
                if kw in typ:
                    candidate = outkey
                    break

        if candidate:
            # try to get residue info from children or the element itself
            # often the element contains sub-elements for residues or partners
            # collect all residue-like children
            found_any = False
            for child in elem.iter():
                if 'residue' in tag_localname(child.tag) or tag_localname(child.tag) in ('partner','atom'):
                    ri = find_residue_info(child)
                    if ri:
                        out[candidate].append(ri)
                        found_any = True
            # fallback: the element itself might represent a single interaction object with residue attributes
            if not found_any:
                ri = find_residue_info(elem)
                if ri:
                    out[candidate].append(ri)

    return out

def main(xml_path, json_path):
    tree = ET.parse(xml_path)
    root = tree.getroot()
    collected = walk_and_collect(root)

    # ensure keys exist
    keys = ['hydrogen_bonds','hydrophobic_contacts','ionic_interactions','water_bridges']
    final = {k: collected.get(k, []) for k in keys}
    # write JSON
    with open(json_path, 'w') as f:
        json.dump(final, f, indent=2)
    print(f"Wrote {json_path} with counts: " + ", ".join(f"{k}={len(final[k])}" for k in keys))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python plip_xml_to_json.py report.xml report.json")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
