#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include <QtGui/QWidget>

namespace Ui {
    class progressBar;
}

class progressBar : public QWidget {
    Q_OBJECT
public:
    progressBar(QWidget *parent = 0);
    progressBar(const int node, const int from, const int to, const int curr, QWidget *parent = 0);
    ~progressBar();


protected:
    void changeEvent(QEvent *e);

private:
    Ui::progressBar *m_ui;
};

#endif // PROGRESSBAR_H
