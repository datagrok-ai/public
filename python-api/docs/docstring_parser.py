from typing import List, Dict, Optional
import re

class DocstringParser:
    """Parser for NumPy-style docstrings.
    
    This parser extracts structured information from NumPy-style docstrings,
    specifically focusing on Parameters, Attributes, and Returns sections.
    """
    
    @staticmethod
    def parse(docstring: str) -> Dict[str, List[Dict[str, str]]]:
        """Parse a NumPy-style docstring and extract structured information.
        
        Parameters
        ----------
        docstring : str
            The docstring to parse
            
        Returns
        -------
        Dict[str, List[Dict[str, str]]]
            Dictionary containing parsed sections with their items.
            Keys are 'parameters', 'attributes', 'returns', and 'docstring'.
            Each item is a dictionary with 'name', 'type', and 'description' keys.
            The 'docstring' key contains the main description text.
        """
        # Initialize result dictionary
        result = {
            'parameters': [],
            'attributes': [],
            'returns': [],
            'docstring': []
        }
        
        # Split docstring into lines and remove empty lines
        lines = [line.strip() for line in docstring.split('\n') if line.strip()]
        
        current_section = None
        current_item = None
        main_description = []
        in_section = False
        
        for line in lines:
            # Check for section headers
            if line in ['Parameters', 'Attributes', 'Returns']:
                # Save any accumulated main description
                if main_description:
                    result['docstring'].append(' '.join(main_description))
                    main_description = []
                
                if current_item:
                    result[current_section].append(current_item)
                current_item = None
                current_section = line.lower()
                in_section = True
                continue
                
            # Skip section separator lines (dashes)
            if line.startswith('---'):
                in_section = True
                continue
                
            # Parse item lines
            if in_section:
                # Try to match parameter/attribute/return item
                match = re.match(r'^([\w\[\], ]+)\s*:\s*([\w\[\]., ]+)(?:\s*,\s*(optional|default=[^,]+))*$', line)
                if match:
                    # If we have a previous item, add it to the result
                    if current_item:
                        result[current_section].append(current_item)
                    
                    # Start new item
                    name = match.group(1).strip()
                    type_info = match.group(2).strip()
                    current_item = {
                        'name': name,
                        'type': type_info,
                        'description': ''
                    }    
                elif current_item and line:
                    # Append to current item's description
                    if current_item['description']:
                        current_item['description'] += ' '
                    current_item['description'] += line
                if not current_item and not match and current_section == 'returns':
                     current_item = {
                        'name': '',
                        'type': line,
                        'description': ''
                    }    
            else:
                # Collect main description text
                main_description.append(line)
        
        # Add the last item if exists
        if current_item:
            result[current_section].append(current_item)
            
        # Add any remaining main description
        if main_description:
            result['docstring'].append(' '.join(main_description))
            
        return result

def parse_docstring(docstring: str) -> Dict[str, List[Dict[str, str]]]:
    """Convenience function to parse a docstring.
    
    Parameters
    ----------
    docstring : str
        The docstring to parse
        
    Returns
    -------
    Dict[str, List[Dict[str, str]]]
        Dictionary containing parsed sections with their items and main description
    """
    return DocstringParser.parse(docstring)

# Example usage:
if __name__ == "__main__":
    example_docstring = """Ensure the model has a valid ID, generating one if necessary.
        
        If the model's ID is None, generates a new UUID v1 and assigns it to the model.
        This is typically used before saving a new model instance to the server.
        
        Returns
        -------
        str
            The model's ID (either existing or newly generated)
        """
    
    result = parse_docstring(example_docstring)
    print("Main description:")
    for desc in result['docstring']:
        print(desc)
        print()
        
    print("Parsed sections:")
    for section, items in result.items():
        if section != 'docstring':  # Skip docstring section as we already printed it
            print(f"\n{section.upper()}:")
            for item in items:
                print(f"  Name: {item['name']}")
                print(f"  Type: {item['type']}")
                print(f"  Description: {item['description']}")
                print() 