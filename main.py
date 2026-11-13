from bakery import assert_equal
from dataclasses import dataclass, field
from drafter import *
from sympy import Matrix
from math import lcm, gcd
from functools import reduce

atomic_masses = {
    "H": 1.008,
    "He": 4.003,
    "Li": 6.94,
    "Be": 9.012,
    "B": 10.81,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "Ne": 20.180,
    "Na": 22.990,
    "Mg": 24.305,
    "Al": 26.982,
    "Si": 28.085,
    "P": 30.974,
    "S": 32.06,
    "Cl": 35.45,
    "Ar": 39.948,
    "K": 39.098,
    "Ca": 40.078,
    "Sc": 44.956,
    "Ti": 47.867,
    "V": 50.942,
    "Cr": 51.996,
    "Mn": 54.938,
    "Fe": 55.845,
    "Co": 58.933,
    "Ni": 58.693,
    "Cu": 63.546,
    "Zn": 65.38,
    "Ga": 69.723,
    "Ge": 72.630,
    "As": 74.922,
    "Se": 78.971,
    "Br": 79.904,
    "Kr": 83.798,
    "Rb": 85.468,
    "Sr": 87.62,
    "Y": 88.906,
    "Zr": 91.224,
    "Nb": 92.906,
    "Mo": 95.95,
    "Tc": 98,
    "Ru": 101.07,
    "Rh": 102.906,
    "Pd": 106.42,
    "Ag": 107.868,
    "Cd": 112.414,
    "In": 114.818,
    "Sn": 118.710,
    "Sb": 121.760,
    "Te": 127.60,
    "I": 126.904,
    "Xe": 131.293,
    "Cs": 132.905,
    "Ba": 137.327,
    "La": 138.905,
    "Ce": 140.116,
    "Pr": 140.908,
    "Nd": 144.242,
    "Pm": 145,
    "Sm": 150.36,
    "Eu": 151.964,
    "Gd": 157.25,
    "Tb": 158.925,
    "Dy": 162.500,
    "Ho": 164.930,
    "Er": 167.259,
    "Tm": 168.934,
    "Yb": 173.045,
    "Lu": 174.967,
    "Hf": 178.49,
    "Ta": 180.948,
    "W": 183.84,
    "Re": 186.207,
    "Os": 190.23,
    "Ir": 192.217,
    "Pt": 195.084,
    "Au": 196.967,
    "Hg": 200.592,
    "Tl": 204.38,
    "Pb": 207.2,
    "Bi": 208.980,
    "Po": 209,
    "At": 210,
    "Rn": 222,
    "Fr": 223,
    "Ra": 226,
    "Ac": 227,
    "Th": 232.038,
    "Pa": 231.036,
    "U": 238.029,
    "Np": 237,
    "Pu": 244,
    "Am": 243,
    "Cm": 247,
    "Bk": 247,
    "Cf": 251,
    "Es": 252,
    "Fm": 257,
    "Md": 258,
    "No": 259,
    "Lr": 266,
    "Rf": 267,
    "Db": 268,
    "Sg": 269,
    "Bh": 270,
    "Hs": 277,
    "Mt": 278,
    "Ds": 281,
    "Rg": 282,
    "Cn": 285,
    "Nh": 286,
    "Fl": 289,
    "Mc": 290,
    "Lv": 293,
    "Ts": 294,
    "Og": 294
}

@dataclass
class Reactant:
    name: str
    mass: float 
    coeff: int = 1


@dataclass
class Product:
    name: str
    coeff: int = 1


@dataclass
class State:
    reactants: list[Reactant] #= field(default_factory=list)
    products: list[Product] #= field(default_factory=list)
    balance: str = ""
    limiting: str = ""
    type: str = ""

@route
def index(state: State) -> Page:
    '''
    This function returns the home page, including the default settings for
    the state dataclass.
    
    Args:
        state (State): represents the basic features of the entire chem
        calculator
    
    Returns:
        Page: Updated page with state values
    '''
    return Page(state, [
        Header("Welcome to Chem Calculator"),

        "Enter a chemical reaction below. Separate "
        "the reactants and products using a '+'. "
        "The reaction does not have to be balanced.",

        "Reactants:",
        TextBox("chem_reactants", "ex. N2+H2"),

        "Products:",
        TextBox("chem_products", "ex. NH3"),

        "Enter the mass of each reactant in grams, "
        "one mass per line, in the same order as "
        "those above:",

        TextArea("masses", "ex. \n28.02\n6.048"),

        Button("Balance Chemical Reaction", "balance_reaction")
    ])


def parse_formula(formula: str):
    '''
    This function uses the dictionary of elements from parse_group
    to return a final dictionary, while confirming the compound
    isn't invalid.
    
    Args:
        formula (str): formula of chemical compound
    
    Returns:
        dictionary: element keys with number of each element next to it
    '''
    def parse_group(index):
        '''
        This function uses the reactants and products entered by
        the user to create a dictionary of elements, auditing
        each appearence.
    
        Args:
        index (int): position in string of compound
    
        Returns:
        dictionary: element keys with number of each element next to it
        '''
        elements = {}
        
        while index < len(formula):
            char = formula[index]

            if char == ')':
                return elements, index + 1

            elif char == '(':
                group_elements, index = parse_group(index + 1)

                start = index
                while index < len(formula) and formula[index].isdigit():
                    index += 1

                multiplier = (
                    int(formula[start:index])
                    if start < index else 1
                )

                for element, count in group_elements.items():
                    elements[element] = (
                        elements.get(element, 0)
                        + count * multiplier
                    )

            elif char.isupper():
                element = char
                index += 1

                if index < len(formula) and formula[index].islower():
                    element += formula[index]
                    index += 1

                start = index
                while index < len(formula) and formula[index].isdigit():
                    index += 1

                count = (
                    int(formula[start:index])
                    if start < index else 1
                )

                elements[element] = elements.get(element, 0) + count

            else:
                raise ValueError(
                    f"Invalid formula character: {char}"
                )

        return elements, index

    result, final_index = parse_group(0)

    if final_index != len(formula):
        raise ValueError("Invalid chemical formula.")

    return result


def count_elements(compounds, coefficients=None):
    '''
    This function calculates the total number of atoms
    of each element across reactants and products, accounting
    for coefficient.
    
    Args:
        compounds (list): list of reactants
        coefficients (list): list of coefficients of compounds
        in a balanced chemical reaction
    
    Returns:
        totals (dictionary): with specific elements matched to
        the total number of the element in the reaction
    '''
    totals = {}

    if coefficients is None:
        coefficients = [1] * len(compounds)

    for compound, coefficient in zip(compounds, coefficients):
        element_counts = parse_formula(compound)

        for element, count in element_counts.items():
            totals[element] = (
                totals.get(element, 0) + count * coefficient
            )

    return totals


def format_side(compounds: list[Reactant], coefficients: list[coeff]):
    '''
    This function takes a list of chemical compounds
    and a corresponding list of their coefficients in a
    balanced equation and outputs a formatted, readable
    chemical equation.
    
    Args:
        compounds (list): list of reactants
        coefficients (list): list of coefficients of compounds
        in a balanced chemical reaction
    
    Returns:
        str: readable chemical reaction
    '''
    result = []

    for compound, coefficient in zip(compounds, coefficients):
        if coefficient == 1:
            result.append(compound)
        else:
            result.append(str(coefficient) + " " + compound)

    return " + ".join(result)


@route
def balance_reaction(
    state: State,
    chem_reactants: str,
    chem_products: str,
    masses: str
) -> Page:
    '''
    This function balances a chemical reaction
    
    Args:
        state (State): general dataclass
        chem_reactants (str): user input of chem_reactants
        chem_products (str): user input of chem_products
        masses (str): masses of reactants from user
        
    Returns:
        Page: new page with new buttons and balanced reaction
    '''
    # Separate the chemical formulas
    reactant_names = [
        x.strip() for x in chem_reactants.split("+")
    ]

    product_names = [
        x.strip() for x in chem_products.split("+")
    ]

    if not all(reactant_names) or not all(product_names):
        raise ValueError("Enter valid reactants and products.")

    if len(set(reactant_names)) != len(reactant_names):
        raise ValueError("Enter each reactant only once.")

    # Read the masses
    mass_values = [
        float(x.strip())
        for x in masses.strip().splitlines()
        if x.strip()
    ]

    if len(mass_values) != len(reactant_names):
        raise ValueError("Enter one mass for each reactant.")

    if any(mass <= 0 for mass in mass_values):
        raise ValueError("All reactant masses must be positive.")

    # Validate formulas
    for compound in reactant_names + product_names:
        parse_formula(compound)

    compounds = reactant_names + product_names

    # Identify all unique elements
    all_elements = set()

    for compound in compounds:
        all_elements.update(parse_formula(compound).keys())

    elements = sorted(all_elements)

    # Create the balancing matrix
    rows = []

    for element in elements:
        row = []

        for compound in reactant_names:
            counts = parse_formula(compound)
            row.append(counts.get(element, 0))

        for compound in product_names:
            counts = parse_formula(compound)
            row.append(-counts.get(element, 0))

        rows.append(row)

    matrix = Matrix(rows)
    solutions = matrix.nullspace()

    if len(solutions) != 1:
        raise ValueError("No unique balancing solution found.")

    coefficients = list(solutions[0])

    # Make coefficients positive
    if all(x < 0 for x in coefficients):
        coefficients = [-x for x in coefficients]

    if any(x <= 0 for x in coefficients):
        raise ValueError("Invalid coefficient solution.")

    # Convert coefficients to whole numbers
    denominator = lcm(*[int(x.q) for x in coefficients])

    coefficients = [
        int(x * denominator) for x in coefficients
    ]

    # Simplify coefficients
    common_factor = reduce(gcd, coefficients)

    coefficients = [
        x // common_factor for x in coefficients
    ]

    # Separate reactant and product coefficients
    reactant_coefficients = coefficients[:len(reactant_names)]
    product_coefficients = coefficients[len(reactant_names):]

    # Verify that the reaction is balanced
    left = count_elements(reactant_names, reactant_coefficients)
    right = count_elements(product_names, product_coefficients)

    if left != right:
        raise ValueError("Reaction balancing verification failed.")

    # Create Reactant and Product objects
    state.reactants = [
        Reactant(name, mass, coefficient)
        for name, mass, coefficient in zip(
            reactant_names, mass_values, reactant_coefficients
        )
    ]

    state.products = [
        Product(name, coefficient)
        for name, coefficient in zip(
            product_names, product_coefficients
        )
    ]

    # Format the balanced equation
    state.balance = (
        format_side(reactant_names, reactant_coefficients)
        + " -> "
        + format_side(product_names, product_coefficients)
    )

    return Page(state, [
        Header("Balanced Chemical Reaction"),

        Text(state.balance),

        Button("Find Limiting Reagent", "limiting_reagent"),

        Button("Predict Reaction Type", "type_of_reaction"),

        Button("Return Home", "return_home")
    ])


# --------------------------------------------------
# 6. CALCULATE MOLAR MASS
# --------------------------------------------------

def calculate_molar_mass(compound: str):
    '''
    This function calculates the molar mass of a compound
    
    Args:
        compound (str): formula of compound
    
    Returns:
        molar_mass(float): molar mass of compound
    '''
    elements = parse_formula(compound)

    molar_mass = 0

    for element, count in elements.items():
        if element not in atomic_masses:
            raise ValueError(
                f"Atomic mass for {element} is not available."
            )

        molar_mass += atomic_masses[element] * count

    return molar_mass


def find_limiting_reagent(reactants: list[Reactant]):
    '''
    This function finds the limiting reagent for a set of reactants
    
    Args:
        reactants (list[Reactant]): list of reactants in chemical reaction
    
    Returns:
        str: name of limiting reagent 
    '''
    ratios = []

    for reactant in reactants:
        moles = reactant.mass / calculate_molar_mass(reactant.name)

        ratios.append(moles / reactant.coeff)

    limiting_index = ratios.index(min(ratios))

    return reactants[limiting_index].name


@route
def limiting_reagent(state: State) -> Page:
    '''
    This function outputs the limiting reagent in the reaction
    
    Args:
        state (State): general program dataclass
    
    Returns:
        Page: updated page with new buttons and user-interaction
    '''
    state.limiting = find_limiting_reagent(state.reactants)

    return Page(state, [
        Header("Limiting Reagent"),

        "Your balanced chemical equation:",
        Text(state.balance),

        "Your limiting reagent is:",
        Text(state.limiting),

        Button(
            "How much more limiting reagent do you need?",
            "find_amount_limiting"
        ),

        Button("Return Home", "return_home")
    ])


def predict_reaction_type(reactants: list[Reactant], products: list[Product]):
    '''
    This function predicts the type of reaction inputted
    
    Args:
        reactants (list[Reactant]): list of reaction's reactants
        products (list[Product]): list of reaction's products
    
    Returns:
        str: type of reaction
    '''
    reactant_names = [r.name for r in reactants]
    product_names = [p.name for p in products]

    # Combustion
    if (
        "O2" in reactant_names
        and "CO2" in product_names
        and "H2O" in product_names
    ):
        return "Combustion"

    # Synthesis
    elif len(reactants) > 1 and len(products) == 1:
        return "Synthesis"

    # Decomposition
    elif len(reactants) == 1 and len(products) > 1:
        return "Decomposition"

    # Single or double replacement
    elif len(reactants) == 2 and len(products) == 2:

        reactant_elemental = any(
            len(parse_formula(r.name)) == 1 for r in reactants
        )

        product_elemental = any(
            len(parse_formula(p.name)) == 1 for p in products
        )

        if reactant_elemental and product_elemental:
            return "Single Replacement"

        elif all(
            len(parse_formula(x)) > 1
            for x in reactant_names + product_names
        ):
            return "Double Replacement"

    return "Unknown or Other Reaction Type"


@route
def type_of_reaction(state: State) -> Page:
    '''
    This function formats the type of reaction to display on the page.
    
    Args:
        state (State): general program dataclass
        
    Returns:
        Page: updated page with updated state values
    '''
    state.type = predict_reaction_type(
        state.reactants,
        state.products
    )

    return Page(state, [
        Header("Reaction Type"),

        "Your balanced chemical equation:",
        Text(state.balance),

        "Predicted reaction type:",
        Text(state.type),

        Button("Learn More About This Reaction", "reaction_type_breakdown"),

        Button("Return Home", "return_home")
    ])


@route
def reaction_type_breakdown(state: State) -> Page:
    '''
    This function outputs an explanation of the type of reaction
    
    Args:
        state (State): general program dataclass
        
    Returns:
        Page: updated page with updated state values
    
    '''
    explanations = {
        "Synthesis":
            "A synthesis reaction occurs when two or more "
            "reactants combine to form one product.",

        "Decomposition":
            "A decomposition reaction occurs when one reactant "
            "breaks down into two or more simpler products.",

        "Single Replacement":
            "A single replacement reaction occurs when one "
            "element replaces another element in a compound.",

        "Double Replacement":
            "A double replacement reaction occurs when two "
            "compounds exchange ions to form two new compounds.",

        "Combustion":
            "A combustion reaction occurs when a substance reacts "
            "with oxygen and releases energy. Hydrocarbon "
            "combustion commonly produces carbon dioxide and water."
    }

    explanation = explanations.get(
        state.type,
        "This reaction does not fit one of the common reaction "
        "types in this program."
    )

    return Page(state, [
        Header("Reaction Explanation"),

        "Your reaction type:",
        Text(state.type),

        "What does this mean?",
        Text(explanation),

        Button("Return Home", "return_home")
    ])

@route
def find_amount_limiting(state: State) -> Page:
    '''
    This function outputs the amount of extra limiting
    reagent needed so the current limiting reagent will not
    be limiting
    
    Args:
        state (State): general program dataclass
        
    Returns:
        Page: updated page with updated state values
        
    '''
    reactants = state.reactants

    if len(reactants) < 2:
        raise ValueError("At least two reactants are required.")

    # Calculate moles of each reactant
    moles = [
        reactant.mass / calculate_molar_mass(reactant.name)
        for reactant in reactants
    ]

    # Calculate mole-to-coefficient ratios
    ratios = [
        mole / reactant.coeff
        for mole, reactant in zip(moles, reactants)
    ]

    # Find the limiting reagent
    limiting_index = ratios.index(min(ratios))

    # Find the smallest ratio among the other reactants
    other_ratios = [
        ratio
        for i, ratio in enumerate(ratios)
        if i != limiting_index
    ]

    target_ratio = min(other_ratios)

    # Calculate how many more moles are needed
    target_moles = (
        reactants[limiting_index].coeff * target_ratio
    )

    additional_moles = max(
        0,
        target_moles - moles[limiting_index]
    )

    # Convert additional moles into grams
    molar_mass = calculate_molar_mass(
        reactants[limiting_index].name
    )

    additional_mass = additional_moles * molar_mass

    return Page(state, [
        Header("Additional Limiting Reagent Needed"),

        "Your limiting reagent is:",
        Text(reactants[limiting_index].name),

        "Additional mass needed to reach the balanced ratio:",
        Text(f"{additional_mass:.2f} grams"),

        Button("Return Home", "return_home")
    ])

@route
def return_home(state: State) -> Page:
    '''
    This function resets the chem calculator by
    going to the home page
    
    Args:
        state (State): general program dataclass
        
    Returns:
        Page: updated home page with new state values
        
    '''
    return index(state)

start_server(State([],[],"","",""))

