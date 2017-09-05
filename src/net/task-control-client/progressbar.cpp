#include "progressbar.h"
#include "ui_progressbar.h"

progressBar::progressBar(QWidget *parent) :    QWidget(parent),  m_ui(new Ui::progressBar) {
    m_ui->setupUi(this);
}

progressBar::progressBar(const int node, const int from, const int to, const int curr, QWidget *parent) :    QWidget(parent),  m_ui(new Ui::progressBar) {
    m_ui->setupUi(this);

    m_ui->bar->setValue(double(curr-from+1)/double(to-from+1)*100);
    QString str;
    m_ui->node->setText("Node: "+str.setNum(node));
    QString f,t;
    m_ui->range->setText(f.setNum(from)+".."+t.setNum(to));
}

progressBar::~progressBar()
{
    delete m_ui;
}

void progressBar::changeEvent(QEvent *e)
{
    QWidget::changeEvent(e);
    switch (e->type()) {
    case QEvent::LanguageChange:
        m_ui->retranslateUi(this);
        break;
    default:
        break;
    }
}
