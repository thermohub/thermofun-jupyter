import ipywidgets as widgets
def log_progress(sequence, every=None, size=None, name='Items', progress_out=widgets.Output()):
    from ipywidgets import IntProgress, HTML, HBox
    from IPython.display import display, clear_output
    global progress
    # global progress_out
    
    with progress_out:
        clear_output(wait=True)

        is_iterator = False
        if size is None:
            try:
                size = len(sequence)
            except TypeError:
                is_iterator = True
        if size is not None:
            if every is None:
                if size <= 200:
                    every = 1
                else:
                    every = int(size / 200)     # every 0.5%
        else:
            assert every is not None, 'sequence is iterator, set every'

        if is_iterator:
            progress = IntProgress(min=0, max=1, value=1)
            progress.layout={'width':'214px', 'max_width':'214px'}
            progress.bar_style = 'info'
        else:
            progress = IntProgress(min=0, max=size, value=0)
            progress.layout={'width':'214px', 'max_width':'214px'}
        label = HTML()
        box = HBox(children=[progress, label])
        display(box)

        index = 0
        try:
            for index, record in enumerate(sequence, 1):
                if index == 1 or index % every == 0:
                    if is_iterator:
                        label.value = '{name}: {index} / ?'.format(
                            name=name,
                            index=index
                        )
                    else:
                        progress.value = index
                        label.value = u'{name}: {index} / {size}'.format(
                            name=name,
                            index=index,
                            size=size
                        )
                yield record
        except:
            progress.bar_style = 'danger'
            raise
        else:
            progress.bar_style = 'success'
            progress.value = index
            if index == 0:
                label.value = "{name}: {index}".format(
                name=name,
                index=str(index or '?')
                )
            else:
                label.value = 'Finished reading database'.format(
                name=name,
                index=str(index or '?')
                )