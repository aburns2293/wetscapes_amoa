import xml.etree.ElementTree as ET
tree = ET.parse('Data/metadata.xml')
root = tree.getroot()
return(root.tag)

