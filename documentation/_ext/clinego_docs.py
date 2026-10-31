"""Local Sphinx directives for the CLINE-GO site.

``in-development``: the one styled block that marks a feature not built yet. Its title and look
live here and in ``_static/clinego.css`` (class ``clinego-indev``), so every page shows it the
same way. In MyST::

    :::{in-development}
    A Shiny form for custom climates is being designed; today a script builds the garden tables.
    :::

``slot``: a dashed placeholder box for content that a later run fills in (a figure, the expected
numbers, a run time). The argument is its title; the body says what will go there::

    :::{slot} Figure: ch3_pca_structure
    PCA of the 336 trees, coloured by sNMF ancestry at the chosen K.
    :::

Wave 2 replaces every slot with the real ``{figure}`` or text; ``grep -rn '{slot}'`` lists what
is left.
"""

from docutils import nodes
from docutils.parsers.rst import Directive
from docutils.parsers.rst.directives.admonitions import BaseAdmonition

IN_DEVELOPMENT_TITLE = "In development"
IN_DEVELOPMENT_CLASS = "clinego-indev"
SLOT_CLASS = "clinego-slot"


class InDevelopment(BaseAdmonition):
    """A generic admonition with a fixed title and the ``clinego-indev`` class."""

    node_class = nodes.admonition
    required_arguments = 0
    optional_arguments = 0
    final_argument_whitespace = False
    has_content = True
    option_spec = {}

    def run(self):
        self.arguments = [IN_DEVELOPMENT_TITLE]
        self.options = {"class": [IN_DEVELOPMENT_CLASS]}
        return super().run()


class Slot(Directive):
    """A titled placeholder box (class ``clinego-slot``)."""

    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True
    has_content = True
    option_spec = {}

    def run(self):
        box = nodes.container(classes=[SLOT_CLASS])
        title = nodes.paragraph(classes=[f"{SLOT_CLASS}-title"])
        title += nodes.strong(text=self.arguments[0])
        box += title
        if self.content:
            self.state.nested_parse(self.content, self.content_offset, box)
        return [box]


def setup(app):
    app.add_directive("in-development", InDevelopment)
    app.add_directive("slot", Slot)
    return {"version": "1", "parallel_read_safe": True, "parallel_write_safe": True}
