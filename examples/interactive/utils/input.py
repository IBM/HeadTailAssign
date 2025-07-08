import ipywidgets as widgets
from ipywidgets import HBox, Label, Layout, IntSlider
from rdkit.Chem.Draw import IPythonConsole
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D, MolsToGridImage
from IPython.display import SVG, display, clear_output
import io
import pandas as pd
import sys

style = {'description_width': 'initial'}
box_layout = Layout(display='inline-flex',
                        flex_flow='row wrap',
                        width='auto')

class Input(widgets.VBox):
    #https://stackoverflow.com/a/61008455

    def __init__(self):
        self.widget_options = widgets.RadioButtons(
        options=['One polymer', 'Multiple polymers'],
        description='How many polymers would you like to run?',
        style=style,
        index=None,
        layout=box_layout,
        disabled=False
    )
        self.bool_widget_holder = widgets.HBox()
        children = [
            self.widget_options,
            self.bool_widget_holder,
        ]

        self.widget_options.observe(self._add_bool_widgets, names=['value'])
        super().__init__(children=children)

        # self.input_value = self._add_bool_widgets

    def _add_bool_widgets(self, widg):
        
        new_widgets = []
        if self.widget_options.value == 'Multiple polymers':
            new_widget = widgets.FileUpload(
                                accept='',  # Accepted file extension e.g. '.txt', '.pdf', 'image/*', 'image/*,.pdf'
                                description='Upload your file',
                                multiple=False  # True to accept multiple files upload else False
    )
            submit_button=widgets.Button(description='View molecules')

            def on_button_clicked(b):
                #https://stackoverflow.com/a/67425117
                input_file = list(self.bool_widget_holder.children[0].value.values())[0]
                content = input_file['content']
                content = io.StringIO(content.decode('utf-8'))
                df = pd.read_csv(content)
                if len(df.index) > 50:
                    display(f"This demo can only run less than 50 monomers per run. Your file has {len(df.index)} monomers")
                    display(f"Please reduce the number of monomers and try again.")
                    sys.exit()
                else:

                    mol = [Chem.MolFromSmiles(m) for m in list(df.monomers_list)]
                    
                    display(MolsToGridImage(mols=mol, legends=list(df.name), subImgSize=(400, 400), molsPerRow=3))
        
            submit_button.on_click(on_button_clicked)
            new_widgets.append(new_widget)
            new_widgets.append(submit_button)
        if self.widget_options.value == 'One polymer':
            new_widget = widgets.Text(
                                    value='C=CC1CCCCC1',
                                    placeholder='Type something',
                                    description='Please, type the SMILES of your monomer:',
                                    style=style,
                                    layout=box_layout,
                                    disabled=False   
    )
            submit_button=widgets.Button(description='Submit text')
            value_list = []

            def on_button_clicked(b):
                # clear_output(wait=True)
                value_list.append(new_widget.value)
                if len(value_list) == 1:
                    mol = Chem.MolFromSmiles(value_list[0])
                    d2d = rdMolDraw2D.MolDraw2DSVG(350,300)
                    d2d.DrawMolecule(mol)
                    d2d.FinishDrawing()

                    display(SVG(data=d2d.GetDrawingText()))
        
            submit_button.on_click(on_button_clicked)
            new_widgets.append(new_widget)
            new_widgets.append(submit_button)
            
        self.bool_widget_holder.children = tuple(new_widgets)

        return new_widget.value

    @property
    def get(self):
        out = self.bool_widget_holder.children
        return out